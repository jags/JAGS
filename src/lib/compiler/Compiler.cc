#include <config.h>

#include <compiler/Compiler.h>
#include <compiler/ParseTree.h>
#include <graph/LogicalNode.h>
#include <graph/ConstantNode.h>
#include <graph/StochasticNode.h>
#include <graph/DevianceNode.h>
#include <graph/AggNode.h>
#include <graph/NodeError.h>
#include <sarray/RangeIterator.h>
#include <util/nainf.h>
#include <util/dim.h>
#include <util/integer.h>

#include "MixCompiler.h"

#include <utility>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <list>
#include <algorithm>
#include <string>
#include <set>

using std::string;
using std::vector;
using std::list;
using std::map;
using std::pair;
using std::invalid_argument;
using std::runtime_error;
using std::logic_error;
using std::ostringstream;
using std::min;
using std::max;
using std::set;
using std::fabs;

#include <sstream>
template<class T>
std::string ToString(const T& val)
{
    ostringstream strm;
    strm << val;
    return strm.str();
}

Node * Compiler::constFromTable(ParseTree const *p)
{
    // Get a constant value directly from the data table

    if (!_index_expression) {
	throw logic_error("Can only call constFromTable inside index expression");
    }

    map<string,SArray>::const_iterator i = _data_table.find(p->name());
    if (i == _data_table.end()) {
	return 0;
    }
    SArray const &sarray = i->second;

    Range subset_range = getRange(p, sarray.range());
    if (isNULL(subset_range)) {
	return 0;
    }
    else {
	// Range expression successfully evaluated
	Node *cnode = 0;
	if (subset_range.length() > 1) {
            
	    RangeIterator i(subset_range);
	    unsigned int n = subset_range.length();
	    //double const *v = sarray.value();
	    vector<double> const &v = sarray.value();
            vector<double> value(n);
	    for (unsigned int j = 0; j < n; ++j, i.nextLeft()) {
		unsigned int offset = sarray.range().leftOffset(i);
		value[j] = v[offset];
		if (value[j] == JAGS_NA) {
                    return 0;
                }
            }
            cnode = new ConstantNode(subset_range.dim(false), value, 
				     _model.nchain());

	}
	else {
	    unsigned int offset = 
		sarray.range().leftOffset(subset_range.lower());  
	    double value = sarray.value()[offset];
	    if (value == JAGS_NA) {
		return 0;
	    }
	    else {
		cnode = new ConstantNode(value, _model.nchain());
	    }
	    return cnode;
	}
	_index_graph.add(cnode);
	return cnode;
    }
}

bool Compiler::indexExpression(ParseTree const *p, int &value)
{
    /* 
       Evaluates an index expression.

       Index expressions occur in three contexts:
       1) In the limits of a "for" loop
       2) On the left hand side of a relation
       3) On the right hand side of a relation
     
       They are scalar, integer-valued, constant expressions.  We
       return true on success and the result is written to the
       parameter value.
    */
    
    /* 
       The counter _index_expression is non-zero if we are inside an
       Index expression. This invokes special rules in the functions
       getParameter and getArraySubset.  The counter tracks the levels
       of nesting of index expressions.
       
       The graph _index_graph holds the Nodes created during the 
       evaluation of the index expression.
    */

    _index_expression++;
    Node *node = getParameter(p);
    _index_expression--;

    Graph tgraph;
    if (_index_expression == 0) {
	if (node) {
	    tgraph.add(node); //Preserves node when _index_graph is cleared
	}
        _index_graph.clear();
    }

    if (!node) {
	return false;
    }
    if (node->length() != 1) {
	throw NodeError(node,"Vector value in index expression"); 
    }
    if (node->isObserved()) {
	bool is_integer = true;
	value = checkInteger(node->value(0)[0], is_integer);
	if (!is_integer) {
	    throw NodeError(node, 
			    "Index expression evaluates to non-integer value");
	}
	return true;
    }
    else {
	return false;
    }
}

Range Compiler::getRange(ParseTree const *p, Range const &default_range)
{
  /* 
     Evaluate a range expression. If successful, it returns the range
     corresponding to the expression.  If unsuccessful (due to missing
     values) returns a null range.
     
     The default_range argument provides default values if the range
     expression is blank: e.g. foo[] or bar[,1].  The default range 
     may be a null range, in which case, missing indices will result in
     failure.
  */
  
    vector<ParseTree*> const &range_list = p->parameters();
    string const &name = p->name();

  if (range_list.empty()) {
      //An empty range expression implies the default range
      return default_range;
  }

  // Check size and integrity of range expression
  unsigned int size = range_list.size();
  if (!isNULL(default_range) && size != default_range.ndim(false)) {
    string msg = string("dimension mismatch taking subset of ") + name;
    throw runtime_error(msg);
  }
  for (unsigned int i = 0; i < size; ++i) {
    if (range_list[i]->treeClass() != P_RANGE) {
      throw logic_error("Malformed parse tree. Expected range expression");
    }
  }
  
  // Now step through and evaluate lower and upper index expressions
  vector<int> lower(size), upper(size);
  for (unsigned int i = 0; i < size; i++) {
    switch (range_list[i]->parameters().size()) {
    case 0:
      // Empty index implies default range
      if (isNULL(default_range)) {
	return default_range;
      }
      lower[i] = default_range.lower()[i];
      upper[i] = default_range.upper()[i];
      break;
    case 1:
      // Single index implies lower == upper
      if (!indexExpression(range_list[i]->parameters()[0], lower[i])) {
	return Range();
      }
      else {
	upper[i] = lower[i];
      }
      break;
    case 2:
      if (!indexExpression(range_list[i]->parameters()[0], lower[i]) ||
	  !indexExpression(range_list[i]->parameters()[1], upper[i])) {
	return Range();
      }
      break;
    default:
      throw logic_error("Malformed parse tree in index expression");
    }
  }
  
  //FIXME: Give informative error message if we request invalid range
  return Range(lower, upper);
}

Range Compiler::VariableSubsetRange(ParseTree const *var)
{
  /*
    Get the range of a subset expression of a variable on the LHS of a
    relation.  This means that the subset expression must be constant.
  */
  if (var->treeClass() != P_VAR) {
    throw logic_error("Expecting variable expression");
  }  

  string const &name = var->name();
  if (_countertab.getCounter(name)) {
    throw runtime_error(string("Counter cannot appear on LHS of relation: ")
			+ name);
  }
  NodeArray *array = _model.symtab().getVariable(name);
  if (array) {
    // It's a declared node
    vector<ParseTree*> const &range_list = var->parameters();
    
    if (range_list.empty()) {
	//Missing range implies the whole node
	return array->range();
    }
    if (range_list.size() != array->range().ndim(false)) {
	throw runtime_error(string("Dimension mismatch in subset expression") +
			    " of variable " + name);
    }
    Range range = getRange(var, array->range());
    if (isNULL(range)) {
	throw runtime_error(string("Missing values in subset expression ") 
			    + "of variable " + name);
    }
	return range;
  }
  else {
      // Undeclared node
      Range range = getRange(var, Range());
      if (isNULL(range)) {
	  throw runtime_error(string("Cannot evaluate subset expression for ")
			      + "undeclared variable " + name);
      }
      return range;
  }
}

Range Compiler::CounterRange(ParseTree const *var)
{
  /* The range expression for a counter differs from that of
     a variable in that it is
     1) one-dimensional
     2) may not be empty
     Further, no variables are created for counters in the
     Symbol Table
  */
  if (var->treeClass() != P_COUNTER) {
    throw logic_error("Expecting counter expression");
  }
  if (var->parameters().size() != 1) {
    throw logic_error("Invalid counter expression");
  }
  Range range();
  
  ParseTree const *prange = var->parameters()[0];
  if (prange->treeClass() != P_RANGE) {
    throw logic_error("Expecting range expression");
  }

  unsigned int size = prange->parameters().size();
  if (size < 1 || size > 2) {
    throw logic_error(string("Invalid range expression for counter")
		      + var->name());
  }
  int lower;
  if(!indexExpression(prange->parameters()[0], lower)) {
    throw runtime_error(string("Unable to evaluate lower index of counter ")
			+ var->name());
  }
  int upper;
  if (prange->parameters().size() == 2) {
    if (!indexExpression(prange->parameters()[1], upper)) {
      throw runtime_error(string("Unable to evaluate upper index of counter ")
			  + var->name());
    }
  }
  else {
    upper = lower;
  }
  
  if (lower > upper) {
    return Range();
  }
  else {
    return Range(vector<int>(1, lower), vector<int>(1, upper));
  }
}

Node * Compiler::getSubsetNode(ParseTree const *var)
{
    if (var->treeClass() != P_VAR) {
	throw logic_error("Expecting variable expression");
    }
    
    NodeArray *array = _model.symtab().getVariable(var->name());
    if (array == 0) {
	if (_strict_resolution) {
	    throw runtime_error(string("Unknown variable ") + var->name());
	}
	else {
	    return 0;
	}
    }
    Range subset_range = getRange(var, array->range());
    if (isNULL(subset_range)) {
	return 0;
    }
    if (!array->range().contains(subset_range)) {
	throw runtime_error(string("Subset ") + var->name() 
			    + print(subset_range)
			    + " out of range");
    }
    Node *node = array->getSubset(subset_range, _model.graph());
    if (node == 0 && _strict_resolution) {
	throw runtime_error(string("Unable to resolve parameter ") + 
			    array->name() + print(subset_range) +
			    " (one of its ancestors may be undefined)");
    }
    return node;
}

Node *Compiler::getArraySubset(ParseTree const *p)
{
    Node *node = 0;
    
    switch(p->treeClass() != P_VAR) {
	throw logic_error("Expecting expression");
    }

    Counter *counter = _countertab.getCounter(p->name()); //A counter
    if (counter) {
	if (_index_expression) {
	    node = new ConstantNode((*counter)[0], _model.nchain());
	    _index_graph.add(node);
	}
	else {
	    node = _constantfactory.getConstantNode((*counter)[0], 
						    _model.graph());
	}
    }
    else {
	NodeArray *array = _model.symtab().getVariable(p->name());
	if (array) {
	    Range subset_range = getRange(p, array->range());
	    if (!isNULL(subset_range)) {
		node = getSubsetNode(p); //A fixed subset
	    }
	    else if (!_index_expression) {
		node = getMixtureNode(p, this); //A stochastic subset
	    } 
	}
	else if (_strict_resolution) {
	    //Give an informative error message in case of failure
	    throw runtime_error(string("Unknown parameter ") + p->name());
	}
	
	if (!node && _index_expression) {
	    //It is possible to evaluate an index expression before
	    //any Nodes are available from the symbol table.
	    node = constFromTable(p);
	}
    }

    return node;
}

static Function const *getLink(ParseTree const *t, FuncTab const &functab)
{
  if (t->treeClass() != P_LINK) {
    throw logic_error("Malformed parse tree: Expected link function");
  }

  Function const *func = functab.findInverse(t->name());
  if (func == 0) {
    string msg("Unable to find inverse of link function ");
    msg.append(t->name());
    throw runtime_error(msg);
  }
  else {
    return func;
  }
}

static Function const *getFunction(ParseTree const *t, FuncTab const &functab)
{
    if (t->treeClass() != P_FUNCTION) 
	throw logic_error("Malformed parse tree: Expected function");

    Function const *func = functab.find(t->name());
    if (func == 0) {
	string msg("Unable to find function ");
	msg.append(t->name());
	throw runtime_error(msg);
    }    

    return func;
}

static Distribution const *getDistribution(ParseTree const *pstoch_rel, 
				           DistTab const &table)
{
  //Get the distribution from a stochastic relation
  
  if (pstoch_rel->treeClass() != P_STOCHREL) {
    throw logic_error("Malformed parse tree. Expecting stochastic relation");
  }
  
  // Get the distribution
  ParseTree const *pdist = pstoch_rel->parameters()[1];
  if (pdist->treeClass() != P_DENSITY) {
    throw logic_error("Malformed parse tree. Expected density expression");
  }
  Distribution const *dist = table.find(pdist->name());
  if (dist == 0) {
    throw runtime_error(string("Unknown distribution: ") + pdist->name());
  }
  return dist;
}

Node *Compiler::getLength(ParseTree const *p, SymTab const &symtab)
{
    if (p->treeClass() != P_LENGTH) {
	throw logic_error("Malformed parse tree. Expecting dim expression");
    }
    ParseTree const *var = p->parameters()[0];
    if (var->treeClass() != P_VAR) {
	throw logic_error("Malformed parse tree. Expecting variable name");
    }
    NodeArray const *array = symtab.getVariable(var->name());
    if (array) {
	Range subset_range = getRange(var, array->range());
	if (isNULL(subset_range)) {
	    return 0;
	}
	else {
	    double length = product(subset_range.dim(true));
	    if (_index_expression) {
		Node *node = new ConstantNode(length, _model.nchain());
		_index_graph.add(node);
	    return node;
	    }
	    else {
		return _constantfactory.getConstantNode(length, _model.graph());
	    }
	}
    }
    else {
	return 0;
    }
}


Node *Compiler::getDim(ParseTree const *p, SymTab const &symtab)
{
    if (p->treeClass() != P_DIM) {
	throw logic_error("Malformed parse tree. Expecting length expression");
    }
    ParseTree const *var = p->parameters()[0];
    if (var->treeClass() != P_VAR) {
	throw logic_error("Malformed parse tree. Expecting variable name");
    }
    NodeArray const *array = symtab.getVariable(var->name());
    if (array) {
	Range subset_range = getRange(var, array->range());
	if (isNULL(subset_range)) {
	    return 0;
	}
	else {
	    vector<unsigned int> idim = subset_range.dim(false);
	    vector<double> ddim(idim.size());
            for (unsigned int j = 0; j < idim.size(); ++j) {
                ddim[j] = idim[j];
	    }

	    vector<unsigned int> d(1, idim.size());

	    if (_index_expression) {
		Node *node = new ConstantNode(d, ddim, _model.nchain());
		_index_graph.add(node);
	    return node;
	    }
	    else {
		return _constantfactory.getConstantNode(d, ddim, 
							_model.graph());
	    }
	}
    }
    else {
	return 0;
    }
}



/*
 * Evaluates the expression t, and returns a pointer to a Node. If the
 * expression cannot be evaluated, a NULL pointer is returned. 
 */
Node * Compiler::getParameter(ParseTree const *t)
{
    vector<Node const *> parents;
    Node *node = 0;

    switch (t->treeClass()) {
    case P_VALUE:
	if (_index_expression) {
	    node = new ConstantNode(t->value(), _model.nchain());
	    _index_graph.add(node);
	}
	else {
	    node =  _constantfactory.getConstantNode(t->value(), 
						     _model.graph());
	}
	break;
    case P_VAR:
	node = getArraySubset(t);
	break;
    case P_LENGTH:
	node = getLength(t,_model.symtab());
	break;
    case P_DIM:
        node = getDim(t, _model.symtab());
	break;
    case P_FUNCTION: case P_LINK:
	if (getParameterVector(t, parents)) {
	    Function const *func = 0;
	    if (t->treeClass() == P_FUNCTION) {
		func = getFunction(t, funcTab());
	    }
	    else {
		func = getLink(t, funcTab());
	    }
	    if (_index_expression) {
		node = new LogicalNode(func, parents);
		_index_graph.add(node);
	    }
	    else {
		node = _logicalfactory.getLogicalNode(func, parents, 
						      _model.graph());
	    }
	}
	break;
    default:
	throw  logic_error("Malformed parse tree.");
	break;
    }

    if (!node)
        return 0;

    /* Initialize deterministic nodes now, if they are functions of 
       data, to aid in node recycling */
    if (!node->isRandomVariable()) {
      node->initializeData();
    }

    if (_index_expression && !node->isObserved()) {
	return 0;
    }
    else {
	return node;
    }
}

/*
 * Before creating the node y <- foo(a,b), or z ~ dfoo(a,b), the parent
 * nodes must a,b be created. This expression evaluates the vector(a,b)
 * Arguments are the same as for getParameter.
 */
bool Compiler::getParameterVector(ParseTree const *t,
				  vector<Node const *> &parents)
{
    if (!parents.empty()) {
	throw logic_error("parent vector must be empty in getParameterVector");
    }

    switch (t->treeClass()) {
    case P_FUNCTION: case P_LINK:
	for (unsigned int i = 0; i < t->parameters().size(); ++i) {
	    Node *node = getParameter(t->parameters()[i]);
	    if (node) {
		parents.push_back(node);
	    }
	    else {
		parents.clear();
		return false;
	    }
			  
	}
	break;
    default:
	throw logic_error("Invalid Parse Tree. Expected function or operator.");
    }
    return true;
}

Node * Compiler::allocateStochastic(ParseTree const *stoch_relation)
{
    ParseTree const *distribution = stoch_relation->parameters()[1];  
    vector<Node const *> parameters;
    Node *lBound = 0, *uBound = 0;

    // Create the parameter vector
    vector<ParseTree*> const &param_list = distribution->parameters();
    for (unsigned int i = 0; i < param_list.size(); ++i) {
	Node *param = getParameter(param_list[i]);
	if (param) {
	    parameters.push_back(param);
	}
	else {
	    return 0;
	}
    }

    // Set upper and lower bounds, if truncated
    if (stoch_relation->parameters().size() == 3) {
	ParseTree const *truncated = stoch_relation->parameters()[2];
	ParseTree const *ll = truncated->parameters()[0];
	ParseTree const *ul = truncated->parameters()[1];
	if (ll) {
	    lBound = getParameter(ll);
	    if (!lBound) {
		return 0;
	    }
	}
	if (ul) {
	    uBound = getParameter(ul);
	    if (!uBound) {
		return 0;
	    }
	}
    }

  
    StochasticNode *snode = 0;
    Distribution const *dist = getDistribution(stoch_relation, distTab());

    //Search data table to see if this is an observed node
    ParseTree *var = stoch_relation->parameters()[0];
    string const &name = var->name();
    map<string,SArray>::const_iterator q = _data_table.find(name);
    if (q != _data_table.end()) {
	vector<double> const &data_value = q->second.value();
        Range const &data_range = q->second.range();

	Range target_range = VariableSubsetRange(var);
	unsigned int nmissing = 0;
	SArray this_data(target_range.dim(true));
	unsigned int i = 0;
	for (RangeIterator p(target_range); !p.atEnd(); p.nextLeft()) {
	    unsigned int j = data_range.leftOffset(p);
	    if (data_value[j] == JAGS_NA) {
		++nmissing;
	    }
	    else {
		this_data.setValue(data_value[j], i++);
	    }
	}

	if (nmissing == 0) {
	    snode = new StochasticNode(dist, parameters, lBound, uBound);
	    if (this_data.dim(true) != drop(snode->dim())) {
		string msg = "Dimension mismatch between node and data:\n";
		msg.append("Expected ");
		msg.append(print(Range(snode->dim())));
		msg.append(" got ");
		msg.append(print(Range(this_data.dim(true))));
		throw NodeError(snode, msg);
	    }
	    snode->setObserved(this_data.value());
	}
	else if (nmissing != this_data.length()) {
	    throw NodeError(snode, "Node cannot be partially missing");
	}
    }
 
    if (snode == 0) {
	snode =  new StochasticNode(dist, parameters, lBound, uBound);
    }

    _model.graph().add(snode);
    return snode;
}

Node * Compiler::allocateLogical(ParseTree const *rel)
{
    ParseTree *expression = rel->parameters()[1];
    Node *node = 0;
    vector <Node const *> parents;

    switch (expression->treeClass()) {
    case P_VALUE: 
	node = new ConstantNode(expression->value(), _model.nchain());
	_model.graph().add(node); 
	/* The reason we aren't using a ConstantFactory here is to ensure
	   that the nodes are correctly named */
	break;
    case P_VAR: case P_FUNCTION: case P_LINK: case P_LENGTH: case P_DIM:
	node = getParameter(expression);
	break;
    default:
	throw logic_error("Malformed parse tree in Compiler::allocateLogical");
    }

    return node;
}

void Compiler::allocate(ParseTree const *rel)
{
    if (_is_resolved[_n_relations])
	return;
    
    Node *node = 0;
    
    switch(rel->treeClass()) {
    case P_STOCHREL:
	node = allocateStochastic(rel);
	break;
    case P_DETRMREL:
	node = allocateLogical(rel);
    break;
    default:
	throw logic_error("Malformed parse tree in Compiler::allocate");
	break;
    }
    
    SymTab &symtab = _model.symtab();
    if (node) {
	ParseTree *var = rel->parameters()[0];
	NodeArray *array = symtab.getVariable(var->name());
	if (!array) {
	    //Undeclared array. It's size is inferred from the dimensions of
	    //the newly created node
	    symtab.addVariable(var->name(), node->dim());
	    array = symtab.getVariable(var->name());
	    array->insert(node, array->range());
	}
	else {
	    // Check if a node is already inserted into this range
	    Range range = VariableSubsetRange(var);
	    if (array->find(range)) {
		throw runtime_error(string("Attempt to redefine node ") +
				    var->name() + print(range));
	    }
	    array->insert(node, range);
	}
	_n_resolved++;
	_is_resolved[_n_relations] = true;
    }
}

void Compiler::setConstantMask(ParseTree const *rel)
{
    ParseTree const *var = rel->parameters()[0];
    string const &name = var->name();
    map<string,vector<bool> >::iterator p = _constant_mask.find(name);
    if (p == _constant_mask.end()) {
	return;
    }
    map<string,SArray>::const_iterator q = _data_table.find(name);
    if (q == _data_table.end()) {
	throw logic_error ("Error in Compiler::setConstantMask");
    }
    Range range = VariableSubsetRange(var);
    Range const &var_range = q->second.range();
    if (!var_range.contains(range)) {
	throw logic_error("Invalid range in Compiler::setConstantMask.");
    }
    vector<bool> &mask = p->second;
    for (RangeIterator i(range); !i.atEnd(); i.nextLeft()) {
	mask[var_range.leftOffset(i)] = false;
    }
}

void Compiler::getArrayDim(ParseTree const *p)
{
  /* 
     Called by traverseTree, this function calculates the size
     of all arrays from the left-hand side of all 
     relations, and stores the results in the map _node_array_ranges.
  */

  ParseTree const *var = p->parameters()[0];
  string const &name = var->name();

  if(var->parameters().empty()) {
      //No index expession => No info on array size
      return;
  }

  Range new_range = VariableSubsetRange(var);

  map<string, vector<vector<int> > >::iterator i = 
      _node_array_ranges.find(name);
  if (i == _node_array_ranges.end()) {
    //Create a new entry
    vector<vector<int> > ivec;
    ivec.push_back(new_range.lower());
    ivec.push_back(new_range.upper());
    _node_array_ranges.insert(pair<const string, vector<vector<int> > >(name,ivec));
  }
  else {
    //Check against the existing entry, and modify if necessary
    unsigned int ndim = i->second[0].size();
    if (new_range.ndim(false) != ndim) {
      throw runtime_error(string("Inconsistent dimensions for array ") + name);
    }
    else {
	for (unsigned int j = 0; j < ndim; ++j) {
	    i->second[0][j] = min(i->second[0][j], new_range.lower()[j]);
	    i->second[1][j] = max(i->second[1][j], new_range.upper()[j]);
      }
    }
  }
}

void Compiler::writeConstantData(ParseTree const *relations)
{
    /* 
       Values supplied in the data table, but which DO NOT
       appear on the left-hand side of a relation, are constants.
       We have to find these values in order to create the 
       constant nodes that form the top level of any graphical
       model.
    */

    //First we set up the constant mask, setting all values to true by
    //default
    map<string, SArray>::const_iterator p;
    for (p = _data_table.begin(); p != _data_table.end(); ++p) {
	pair<string, vector<bool> > apair;
	apair.first = p->first;
	apair.second = vector<bool>(p->second.length(), true);
	_constant_mask.insert(apair);
    }

    //Now traverse the parse tree, setting node array subsets that
    //correspond to the left-hand side of any relation to be false
    traverseTree(relations, &Compiler::setConstantMask);

    //Create a temporary copy of the data table containing only
    //data for constant nodes
    map<string, SArray> temp_data_table = _data_table;
    map<string, SArray>::iterator p2;
    for(p2 = temp_data_table.begin(); p2 != temp_data_table.end(); ++p2) {
	string const &name = p2->first;
	SArray &temp_data = p2->second;
	vector<bool> const &mask = _constant_mask.find(name)->second;
	for (unsigned long i = 0; i < temp_data.length(); ++i) {
	    if (!mask[i]) {
		temp_data.setValue(JAGS_NA, i);
	    }
	}
    }

    _model.symtab().writeData(temp_data_table);
}

void Compiler::writeRelations(ParseTree const *relations)
{
  writeConstantData(relations);

  // Set up boolean vector for nodes to indicate whether they are
  // resolved or not.
  _is_resolved = new bool[_n_relations];
  for (unsigned int i = 0; i < _n_relations; ++i) {
    _is_resolved[i] = false;
  }
  
  for (unsigned long N = _n_relations; N > 0; N -= _n_resolved) {
    _n_resolved = 0;
    traverseTree(relations, &Compiler::allocate);
    if (_n_resolved == 0) {
      // Try again, but this time throw an exception from getSubsetNode
      _strict_resolution = true;
      traverseTree(relations, &Compiler::allocate);
      // If that didn't work (but it should!) just throw a generic message
      throw runtime_error("Unable to resolve relations");
    }
  }
  delete [] _is_resolved; _is_resolved = 0;
}

void Compiler::traverseTree(ParseTree const *relations, CompilerMemFn fun,
			    bool resetcounter)
{
  /* 
     Traverse parse tree, expanding FOR loops and applying function
     fun to relations.
  */

  if (resetcounter) {
    _n_relations = 0;
  }

  vector<ParseTree*> const &relation_list = relations->parameters();
  for (vector<ParseTree*>::const_iterator p = relation_list.begin(); 
       p != relation_list.end(); ++p) 
    {  
      Counter *counter;
      ParseTree *var;
      
      switch ((*p)->treeClass()) {
      case P_FOR:
	var = (*p)->parameters()[0];
	if (!isNULL(CounterRange(var))) {
	  counter = _countertab.pushCounter(var->name(), CounterRange(var));
	  for (; !counter->atEnd(); counter->next()) {
	    traverseTree((*p)->parameters()[1], fun, false);
	  }
	  _countertab.popCounter();
	}
	break;
      case P_STOCHREL: case P_DETRMREL:
	(this->*fun)(*p);
	_n_relations++;
	break;
      default:
	throw logic_error("Malformed parse tree in Compiler::traverseTree");
	break;
      }
    }
}

Compiler::Compiler(BUGSModel &model, map<string, SArray> const &data_table)
    : _model(model), _countertab(), 
      _data_table(data_table), _n_resolved(0), 
      _n_relations(0), _is_resolved(0), _strict_resolution(false),
      _index_expression(0), _index_graph(),
      _constantfactory(model.nchain())
{
  if (_model.graph().size() != 0)
    throw invalid_argument("Non empty graph in Compiler constructor");
  if (_model.symtab().size() != 0)
    throw invalid_argument("Non empty symtab in Compiler constructor");
}

void Compiler::declareVariables(vector<ParseTree*> const &dec_list)
{
  vector<ParseTree*>::const_iterator p;
  for (p = dec_list.begin() ; p != dec_list.end(); ++p) {
    if ((*p)->treeClass() != P_VAR) {
      throw invalid_argument("Expected variable expression");
    }
  }
  
  for (p = dec_list.begin() ; p != dec_list.end(); ++p) {
    ParseTree const *node_dec = *p;
    string const &name = node_dec->name();
    unsigned int ndim = node_dec->parameters().size();
    if (ndim == 0) {
	// Variable is scalar
	_model.symtab().addVariable(name, vector<unsigned int>(1,1));
    }
    else {
      // Variable is an array
	vector<unsigned int> dim(ndim);
	for (unsigned int i = 0; i < ndim; ++i) {
	    int dim_i;
	    if (!indexExpression(node_dec->parameters()[i], dim_i)) {
		throw runtime_error(string("Unable to calculate dimensions of node ")
				    + name);
	    }
	    if (dim_i <= 0) {
		throw runtime_error(string("Non-positive dimension for node") 
				    + name);
	    }
	    dim[i] = static_cast<unsigned int>(dim_i);
	}
	_model.symtab().addVariable(name, dim);
    }
  }
}

void Compiler::undeclaredVariables(ParseTree const *prelations)
{
    // Get undeclared variables from data table
    map<string, SArray>::const_iterator p = _data_table.begin();
    for (; p != _data_table.end(); ++p) {
	string const &name = p->first;
	NodeArray const *array = _model.symtab().getVariable(name);
	if (array) {
	    if (p->second.range() != array->range()) {
		string msg = string("Dimensions of ") + name + 
		    " in declaration (" + print(array->range()) + 
		    ") conflict with dimensions in data (" + 
		    print(p->second.range()) + ")";
		throw runtime_error(msg);
	    }
	}
	else {
	    _model.symtab().addVariable(name, p->second.dim(false));
	}
    }

  
    // Infer the dimension of remaining nodes from the relations
    traverseTree(prelations, &Compiler::getArrayDim);
    map<string, vector<vector<int> > >::const_iterator i = 
	_node_array_ranges.begin(); 
    for (; i != _node_array_ranges.end(); ++i) {
	if (_model.symtab().getVariable(i->first)) {
	    //Node already declared. Check consistency 
	    NodeArray const * array = _model.symtab().getVariable(i->first);
	    vector<int> const &upper = array->range().upper();
	    if (upper.size() != i->second[1].size()) {
		string msg = "Dimension mismatch between data and model for node ";
		msg.append(i->first);
		throw runtime_error(msg);
	    }
	    for (unsigned int j = 0; j < upper.size(); ++j) {
		if (i->second[1][j] > upper[j]) {
		    string msg =  string("Index out of range for node ") + i->first;
		    throw runtime_error(msg);
		}
	    } 
	}
	else {
	    //Node not declared. Use inferred size
	    vector<int> const &upper = i->second[1];
	    unsigned int ndim = upper.size();
	    vector<unsigned int> dim(ndim);
	    for (unsigned int j = 0; j < ndim; ++j) {
		if (upper[j] <= 0) {
		    string msg = string("Invalid index for node ") + i->first;
		    throw runtime_error(msg);
		}
		else {
		    dim[j] = static_cast<unsigned int>(upper[j]);
		}
	    }
				 
	    _model.symtab().addVariable(i->first, dim);
	}
    }
}

DistTab &Compiler::distTab()
{
    static DistTab _disttab;
    return _disttab;
}

FuncTab &Compiler::funcTab()
{
    static FuncTab _functab;
    return _functab;
}

MixtureFactory& Compiler::mixtureFactory()
{
   return _mixfactory;
}

BUGSModel &Compiler::model() const
{
    return _model;
}
