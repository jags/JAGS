#include <config.h>

#include <compiler/Compiler.h>
#include <compiler/ParseTree.h>
#include <graph/ScalarLogicalNode.h>
#include <graph/VectorLogicalNode.h>
#include <graph/ArrayLogicalNode.h>
#include <graph/LinkNode.h>
#include <graph/ConstantNode.h>
#include <graph/ScalarStochasticNode.h>
#include <graph/VectorStochasticNode.h>
#include <graph/ArrayStochasticNode.h>
#include <graph/AggNode.h>
#include <graph/NodeError.h>
#include <sarray/SimpleRange.h>
#include <sarray/RangeIterator.h>
#include <function/FunctionPtr.h>
#include <distribution/DistPtr.h>
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
#include <sstream>
#include <map>

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
using std::max_element;



namespace jags {

    typedef pair<vector<unsigned long>, vector<double> > CNodeKey;

    static bool lt(CNodeKey const &k1, CNodeKey const &k2) {
	//Sort first on dimension
	if (k1.first < k2.first) {
	    return true;
	}
	else if (k1.first > k2.first) {
	    return false;
	}
	else {
	    //Fuzzy sort on values
	    return lt(&k1.second[0], &k2.second[0], k1.second.size());
	}

    }
    typedef map<CNodeKey, ConstantNode *, fuzzy_less<CNodeKey> > CNodeMap;

#include <sstream>
template<class T> 
string ToString(const T& val)
{
    ostringstream strm;
    strm << val;
    return strm.str();
}

static void CompileError(ParseTree const *p, string const &msg1, 
			 string const &msg2 = "")
{
    string msg = string("Compilation error on line ") + ToString(p->line())
	+ ".";
    if (!msg1.empty()) {
	msg.append("\n");
	msg.append(msg1);
    }
    if (!msg2.empty()) {
	msg.append(" ");
	msg.append(msg2);
    }
    throw runtime_error(msg);
}

void Compiler::getLHSVars(ParseTree const *relations)
{
    /* 
       Utility function to get the names of variables used on the
       left hand side of a relation. These are added to the set
       _lhs_vars. The main use of this information is to give better
       diagnostic messages when we cannot resolve a variable name.

       Although this function has a similar recursive structure to
       traverseTree it expands for loops without defining a counter.
    */
    
    vector<ParseTree*> const &relation_list = relations->parameters();
    for (vector<ParseTree*>::const_reverse_iterator p = relation_list.rbegin(); 
	 p != relation_list.rend(); ++p) 
    {
	TreeClass tc = (*p)->treeClass();
	if (tc == P_STOCHREL || tc == P_DETRMREL) {
	    ParseTree *var = (*p)->parameters()[0];
	    _lhs_vars.insert(var->name());
	}
	else if (tc == P_FOR) {
	    getLHSVars((*p)->parameters()[1]);
	}
    }
}


    
Node * Compiler::constFromTable(ParseTree const *p)
{
    // Get a constant value directly from the data table

    if (!_index_expression) {
	throw logic_error("Can only call constFromTable inside index expression");
    }

    if (_countertab.getCounter(p->name())) {
	// DO NOT get value from the data table if there is a counter
	// that shares the same name.
	return 0;
    }
    map<string,SArray>::const_iterator i = _data_table.find(p->name());
    if (i == _data_table.end()) {
	return 0;
    }
    SArray const &sarray = i->second;

    Range subset_range = getRange(p, sarray.range());
    if (isNULL(subset_range)) {
	//Range expression not evaluated
	return 0;
    }
    else if (subset_range.length() > 1) {
	//Multivariate constant
	RangeIterator r(subset_range);
	unsigned long n = subset_range.length();
	//double const *v = sarray.value();
	vector<double> const &v = sarray.value();
	vector<double> value(n);
	for (unsigned int j = 0; j < n; ++j, r.nextLeft()) {
	    unsigned long offset = sarray.range().leftOffset(r);
	    value[j] = v[offset];
	    if (value[j] == JAGS_NA) {
		return 0;
	    }
	}
	return getConstant(subset_range.dim(false), value, 
			   _model.nchain(), true);
	
    }
    else {
	//Scalar constant
	unsigned long offset = sarray.range().leftOffset(subset_range.first());  
	double value = sarray.value()[offset];
	if (value == JAGS_NA) {
	    return 0;
	}
	else {
	    return getConstant(value, _model.nchain(), true);
	}
    }
}

bool Compiler::indexExpression(ParseTree const *p, vector<unsigned long> &value)
{
    /* 
       Evaluates an index expression.

       Index expressions occur in three contexts:
       1) In the limits of a "for" loop
       2) On the left hand side of a relation
       3) On the right hand side of a relation
     
       They are integer-valued, constant expressions.  We return true
       on success and the result is written to the parameter value.
    */
    
    /* 
       The counter _index_expression is non-zero if we are inside an
       Index expression. This invokes special rules in the functions
       getParameter and getArraySubset.  The counter tracks the levels
       of nesting of index expressions.
       
       The vector _index_nodes holds the Nodes created during the 
       evaluation of the index expression.
    */

    bool isfixed = true;

    _index_expression++;
    Node *node = getParameter(p);
    _index_expression--;

    if (!node || !node->isFixed()) {
	isfixed = false;
    }
    else {
	for (unsigned int i = 0; i < node->length(); ++i) {
	    double v = node->value(0)[i];
	    if (!checkULong(v)) {
		throw NodeError(node, 
				"Index expression evaluates to invalid value");
	    }
	    value.push_back(asULong(v));
	}
    }

    /* Make sure we always clean up at the top-level call to IndexExpression */

    if (_index_expression == 0) {
	while(!_index_nodes.empty()) {
	    Node *inode = _index_nodes.back();
	    _index_nodes.pop_back();
	    inode->unlinkParents();
	    delete inode;
	}
    }

    return isfixed;
}

Range Compiler::getRange(ParseTree const *p, 
			 SimpleRange const &default_range)
{
  /* 
     Evaluate a range expression. If successful, it returns the range
     corresponding to the expression.  If unsuccessful (due to missing
     values) returns a null range.
     
     The default_range argument provides default values if the range
     expression is blank: e.g. foo[] or bar[,1].  The default range 
     may be a NULL range, in which case, missing indices will result in
     failure.
  */
  
    vector<ParseTree*> const &range_list = p->parameters();
    string const &name = p->name();

    if (range_list.empty()) {
	//An empty range expression implies the default range
	return default_range;
    }

  // Check size and integrity of range expression
  unsigned long size = range_list.size();
  if (!isNULL(default_range) && size != default_range.ndim(false)) {
      CompileError(p, "Dimension mismatch taking subset of", name);
  }
  for (unsigned int i = 0; i < size; ++i) {
    if (range_list[i]->treeClass() != P_RANGE) {
      throw logic_error("Malformed parse tree. Expected range expression");
    }
  }
  
  // Now step through and evaluate lower and upper index expressions
  vector<vector<unsigned long> > scope(size);
  for (unsigned int i = 0; i < size; i++) {
    switch (range_list[i]->parameters().size()) {
    case 0:
      // Empty index implies default range
      if (isNULL(default_range)) {
	  return Range();
      }
      scope[i] = default_range.scope()[i];
      break;
    case 1:
      if (!indexExpression(range_list[i]->parameters()[0], scope[i])) {
	  return Range();
      }
      if (scope[i].empty()) {
	  CompileError(p, "Invalid range");
      }
      break;
    default:
      throw logic_error("Malformed parse tree in index expression");
    }
  }

  if (!isNULL(default_range)) {
      // If a default range is given, the subset cannot be outside of it
      if (!default_range.contains(Range(scope))) {
	  CompileError(p, "Index out of range taking subset of ", name);
      }
  }
  
  return Range(scope);
}

SimpleRange Compiler::VariableSubsetRange(ParseTree const *var)
{
    /*
      Get the range of a subset expression of a variable on the LHS of a
      relation.  This means that the subset expression must be constant
      and it must be a simple range.
    */
    if (var->treeClass() != P_VAR) {
	throw logic_error("Expecting variable expression");
    }  

    string const &name = var->name();
    if (_countertab.getCounter(name)) {
	CompileError(var, "Attempt to redefine counter defined in for loop:",
		     name);
    }
    NodeArray *array = _model.symtab().getVariable(name);
    SimpleRange default_range;
    if (array) {
	// Declared node
	vector<ParseTree*> const &range_list = var->parameters();
    
	if (range_list.empty()) {
	    //Missing range implies the whole node
	    return array->range();
	}
	if (range_list.size() != array->range().ndim(false)) {
	    CompileError(var, "Dimension mismatch in subset expression of",
			 name);
	}
	default_range = array->range();
    }

    Range range = getRange(var, default_range);
    if (isNULL(range)) {
	CompileError(var, "Cannot evaluate subset expression for", name);
    }
    
    //New in 4.1.0: Enforce use of simple ranges on the LHS of a relation
    for (unsigned int i = 0; i < range.ndim(false); ++i) {
	vector<unsigned long> const &indices = range.scope()[i];
	for (unsigned long j = 1; j < indices.size(); ++j) {
	    if (indices[j] != indices[j-1] + 1) {
		string msg = string("Invalid subset expression for ") + name +
		    "\nIndex expressions on the left hand side of a relation"
		    + "\nmust define a contiguous, increasing set of indices";
		CompileError(var, msg);
	    }
	}
    }
    return SimpleRange(range.first(), range.last());
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
  Range range;
  
  ParseTree const *prange = var->parameters()[0];
  if (prange->treeClass() != P_RANGE) {
    throw logic_error("Expecting range expression");
  }

  unsigned long size = prange->parameters().size();
  if (size != 1) {
    throw logic_error(string("Invalid range expression for counter ")
		      + var->name());
  }
  vector<unsigned long> indices;
  if(!indexExpression(prange->parameters()[0], indices)) {
      CompileError(var, "Cannot evaluate range of counter", var->name());
  }
  if (indices.empty()) {
    return Range();
  }
  else {
      return Range(vector<vector<unsigned long> >(1, indices));
  }
}

    Node * 
    Compiler::getConstant(double value, unsigned int nchain, bool observed) 
    {
	return getConstant(vector<unsigned long>(1, 1), vector<double>(1, value),
			   nchain, observed);
    }

    Node * 
    Compiler::getConstant(vector<unsigned long> const &dim, 
			  vector<double> const &value,
			  unsigned int nchain, bool observed)
    {
	ConstantNode * cnode = 0;
	if (_index_expression) {
	    cnode = new ConstantNode(dim, value, nchain, observed);
	    _index_nodes.push_back(cnode);
	}
	else {
	    CNodeMap::key_type k(dim, value);
	    CNodeMap::const_iterator p = _cnode_map.find(k);
	    if (p != _cnode_map.end()) {
		cnode = p->second;
	    }
	    else {
		cnode = new ConstantNode(dim, value, nchain, observed);
		_cnode_map[k] = cnode;
		_model.addNode(cnode);
	    }
	}
	return cnode;
    }

Node *Compiler::getArraySubset(ParseTree const *p)
{
    Node *node = 0;
    
    if (p->treeClass() != P_VAR) {
	throw logic_error("Expecting expression");
    }

    if (_index_expression) {
	//It is possible to evaluate an index expression before
	//any Nodes are available from the symbol table.
	node = constFromTable(p);
	if (node) return node;
    }
    
    Counter *counter = _countertab.getCounter(p->name()); //A counter
    if (counter) {
	node = getConstant((*counter)[0], _model.nchain(), false);
    }
    else {
	NodeArray *array = _model.symtab().getVariable(p->name());
	if (array) {
	    Range subset_range = getRange(p, array->range());
	    if (!isNULL(subset_range)) {
		//A fixed subset
		if (!array->range().contains(subset_range)) {
		    CompileError(p, "Subset out of range:", array->name() +
				 printRange(subset_range));
		}
		node = array->getSubset(subset_range, _model);
		if (node == 0 && _compiler_mode == COLLECT_UNRESOLVED) {
		    /* Nake a note of all subsets that could not be
		       resolved.

		       Example: we have just failed to resolve x[1:5]
		       but we need to know which of x[1], x[2], x[3],
		       x[4], x[5] are empty.
		    */
		    bool empty = true;
		    for (RangeIterator r(subset_range); !r.atEnd();
			 r.nextLeft())
		    {
			if (array->getSubset(SimpleRange(r,r), _model)) {
			    empty = false;
			    break;
			}
		    }
		    if (empty) {
			//No elements found, so report the whole range
			//as missing, e.g. x[1:5]
			pair<string, Range> upair(p->name(), subset_range);
			if (_umap.find(upair) == _umap.end()) {
			    set<unsigned long> lines;
			    lines.insert(p->line());
			    _umap[upair] = lines;
			}
			else {
			    _umap.find(upair)->second.insert(p->line());
			}
		    }
		    else {
			//Some elements found, so report only the missing
			//elements, e.g. x[5]
			for (RangeIterator r(subset_range); !r.atEnd();
			     r.nextLeft())
			{
			    SimpleRange sr(r,r);
			    if (!array->getSubset(sr, _model)) {
				pair<string, Range> upair(p->name(), sr);
				if (_umap.find(upair) == _umap.end()) {
				    set<unsigned long> lines;
				    lines.insert(p->line());
				    _umap[upair] = lines;
				}
				else {
				    _umap.find(upair)->second.insert(p->line());
				}
			    }
			}
		    }
		}
	    }
	    else if (!_index_expression) {
		//A stochastic subset
		node = getMixtureNode(p, this);
		if (node == 0 && _compiler_mode == COLLECT_UNRESOLVED) {
		    getMissingMixParams(p, _umap, this);
		}
	    } 
	}
	else if (_lhs_vars.find(p->name()) == _lhs_vars.end()) {
	    string msg = string("Unknown variable ") + p->name() + "\n" +
		"Either supply values for this variable with the data\n" +
		"or define it  on the left hand side of a relation.";
	    CompileError(p, msg);
	}
	else if (_compiler_mode == GET_DIMS) {
	    string msg = string("Variable `") + p->name() +
		"` cannot be used in an index expression\n" +
		"Try defining it in a data block or supplying data values.";
	}
	else if (_compiler_mode == COLLECT_UNRESOLVED) {
	    string msg = string("Possible directed cycle involving variable `")
		+ p->name() + "` ";
	    CompileError(p, msg);
	}
    }
    return node;
}

static FunctionPtr const &
getFunction(ParseTree const *t, FuncTab const &functab)
{
    if (t->treeClass() != P_FUNCTION) 
	throw logic_error("Malformed parse tree: Expected function");

    FunctionPtr const &func = functab.find(t->name());
    if (isNULL(func)) {
	CompileError(t, "Unknown function:", t->name());
    }    

    return func;
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
	    return getConstant(length, _model.nchain(), false);
	}
    }
    else {
	return 0;
    }
}


Node *Compiler::getDim(ParseTree const *p, SymTab const &symtab)
{
    if (p->treeClass() != P_DIM) {
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
	    vector<unsigned long> idim = subset_range.dim(false);
	    vector<double> ddim(idim.size());
            for (unsigned long j = 0; j < idim.size(); ++j) {
                ddim[j] = idim[j];
	    }

	    vector<unsigned long> d(1, idim.size());
	    return getConstant(d, ddim, _model.nchain(), false);
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
	node = getConstant(t->value(), _model.nchain(), false);
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
    case P_LINK:
	if (getParameterVector(t, parents)) {
	    LinkFunction const *link = funcTab().findLink(t->name());
	    if (!link) {
		CompileError(t, "Unknown link function:", t->name());
	    }
	    node = _logicalfactory.getNode(FunctionPtr(link), parents, _model);
	}
	break;
    case P_FUNCTION:
	if (getParameterVector(t, parents)) {
	    FunctionPtr const &func = getFunction(t, funcTab());
	    if (_index_expression) {
		node = LogicalFactory::newNode(func, parents, _model.nchain());
		_index_nodes.push_back(node);
	    }
	    else {	    
		node = _logicalfactory.getNode(func, parents, _model);
	    }
	}
	break;
    case P_RANGE: case P_BOUNDS:  case P_COUNTER: case P_DENSITY:
    case P_STOCHREL: case P_DETRMREL: case P_FOR: case P_RELATIONS:
    case P_VECTOR: case P_ARRAY: case P_SUBSET: case P_INTERVAL:
	throw  logic_error("Malformed parse tree.");
    }

    if (!node)
        return 0;

    if (_index_expression) {
	//Random variables in index expressions must be observed
	if (node->randomVariableStatus() == RV_TRUE_UNOBSERVED)
	    return 0;
    }

    return node;
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

    bool ok = true;
    switch (t->treeClass()) {
    case P_FUNCTION: case P_LINK: case P_DENSITY:
	if (t->parameters().size() == 0)
	    CompileError(t, "Parameter(s) missing for", t->name());
	for (unsigned int i = 0; i < t->parameters().size(); ++i) {
	    Node *node = getParameter(t->parameters()[i]);
	    if (node) {
		parents.push_back(node);
	    }
	    else {
		ok = false;
	    }
	}
	if (!ok) {
	    parents.clear();
	    return false;
	}
	break;
    case P_VAR: case P_RANGE: case P_BOUNDS:  case P_COUNTER: case P_VALUE:
    case P_STOCHREL: case P_DETRMREL: case P_FOR: case P_RELATIONS:
    case P_VECTOR: case P_ARRAY: case P_DIM: case P_LENGTH: case P_SUBSET:
    case P_INTERVAL:
	throw logic_error("Invalid Parse Tree.");
    }
    return true;
}

Node * Compiler::allocateStochastic(ParseTree const *stoch_relation)
{
    ParseTree const *distribution = stoch_relation->parameters()[1];  

    // Create the parameter vector
    vector<Node const *> parameters;
    if (!getParameterVector(distribution, parameters)) {
	return 0;
    }

    // Set upper and lower bounds
    Node *lBound = 0, *uBound = 0;
    if (stoch_relation->parameters().size() == 3) {
	//Truncated distribution
	ParseTree const *truncated = stoch_relation->parameters()[2];
	if (truncated->treeClass() != P_BOUNDS &&
	    truncated->treeClass() != P_INTERVAL)
	{
	    throw logic_error("Invalid parse tree");
	}
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

    /* 
       Check data table to see if this is an observed node.  If it is,
       we put the data in an array of doubles pointed to by this_data,
       and set data_length equal to the length of the array
    */
    double *this_data = 0;
    unsigned long data_length = 0;

    ParseTree *var = stoch_relation->parameters()[0];
    map<string,SArray>::const_iterator q = _data_table.find(var->name());
    if (q != _data_table.end()) {

	vector<double> const &data_value = q->second.value();
	SimpleRange const &data_range = q->second.range();

	SimpleRange target_range = VariableSubsetRange(var);
	data_length = target_range.length();
	this_data = new double[data_length];

	unsigned int i = 0;
	unsigned int nmissing = 0;
	for (RangeIterator p(target_range); !p.atEnd(); p.nextLeft()) {
	    unsigned long j = data_range.leftOffset(p);
	    if (data_value[j] == JAGS_NA) {
		++nmissing;
	    }
	    this_data[i++] = data_value[j];
	}
	if (nmissing == data_length) {
	    delete [] this_data;
	    this_data = 0;
	    data_length = 0;
	}
	else if (nmissing != 0) {
	    delete [] this_data;
	    CompileError(var, var->name() + printRange(target_range),
			 "is partly observed and partly missing");
	}
    }

    // Check that distribution exists
    string const &distname = distribution->name();
    DistPtr const &dist = distTab().find(distname);
    if (isNULL(dist)) {
	CompileError(distribution, "Unknown distribution:", distname);
    }

    if (!this_data) {
	/* 
	   Special rule for observable functions, which exist both as
	   a Function and a Distribution.  If the node is unobserved,
	   and we find a function matched to the distribution in
	   obsFuncTab, then we create a Logical Node instead.
	*/
	FunctionPtr const &func = obsFuncTab().find(dist);
	if (!isNULL(func)) {
	    LogicalNode *lnode = LogicalFactory::newNode(func, parameters,
							 _model.nchain());
	    _model.addNode(lnode);
            return lnode;
	}
    }	

    /* 
       We allow BUGS-style interval censoring notation for
       compatibility but only allow it if there are no free parameters
       in the distribution
    */
    if (stoch_relation->parameters().size() == 3) {
	ParseTree const *t = stoch_relation->parameters()[2];
	if (t->treeClass() == P_INTERVAL) {
	    for (unsigned int i = 0; i < parameters.size(); ++i) {
		if (!parameters[i]->isFixed()) {
		    CompileError(stoch_relation,
				 "BUGS I(,) notation is only allowed if",
				 "all parameters are fixed");
		}
	    }
	}
    }
    StochasticNode *snode = 0;
    if (SCALAR(dist)) {
	snode =  new ScalarStochasticNode(SCALAR(dist), _model.nchain(),
					  parameters, 
					  lBound, uBound);
    }
    else if (VECTOR(dist)) {
	snode = new VectorStochasticNode(VECTOR(dist), _model.nchain(),
					 parameters,
					 lBound, uBound);
    }
    else if (ARRAY(dist)) {
	snode = new ArrayStochasticNode(ARRAY(dist), _model.nchain(),
					parameters,
					lBound, uBound);
    }
    else {
	throw logic_error("Unable to classify distribution");
    }
    _model.addNode(snode);
    
    if (this_data) {
	snode->setData(this_data, data_length);
	delete [] this_data;
    }
    
    return snode;
}

Node * Compiler::allocateLogical(ParseTree const *rel)
{
    ParseTree *expression = rel->parameters()[1];
    Node *node = 0;
    vector <Node const *> parents;

    switch (expression->treeClass()) {
    case P_VALUE: 
	node = getConstant(expression->value(), _model.nchain(), false);
	break;
    case P_VAR: case P_FUNCTION: case P_LINK: case P_LENGTH: case P_DIM:
	node = getParameter(expression);
	break;
    case P_RANGE: case P_BOUNDS: case P_DENSITY: case P_COUNTER:
    case P_STOCHREL: case P_DETRMREL: case P_FOR: case P_RELATIONS:
    case P_VECTOR: case P_ARRAY: case P_SUBSET: case P_INTERVAL:
	throw logic_error("Malformed parse tree in Compiler::allocateLogical");
    }

    /*
      Check that there are no values in the data table corresponding to
      this node.
    */
    ParseTree *var = rel->parameters()[0];
    map<string,SArray>::const_iterator q = _data_table.find(var->name());
    if (q != _data_table.end()) {
	vector<double> const &data_value = q->second.value();
	SimpleRange const &data_range = q->second.range();
	SimpleRange target_range = VariableSubsetRange(var);

	for (RangeIterator p(target_range); !p.atEnd(); p.nextLeft()) {
	    unsigned long j = data_range.leftOffset(p);
	    if (data_value[j] != JAGS_NA) {
		CompileError(var, var->name() + printRange(target_range),
			     "is a logical node and cannot be observed");
	    }
	}
    }

    return node;
}

void Compiler::allocate(ParseTree const *rel)
{
    if (_is_resolved[_n_relations])
	return;
    
    Node *node = 0;
    
    if (rel->treeClass() == P_STOCHREL) {
	node = allocateStochastic(rel);
    }
    else if (rel->treeClass() == P_DETRMREL) {
	node = allocateLogical(rel);
    }
    else {
	throw logic_error("Malformed parse tree in Compiler::allocate");
    }
    
    SymTab &symtab = _model.symtab();
    if (node) {
	ParseTree *var = rel->parameters()[0];
	NodeArray *array = symtab.getVariable(var->name());
	if (!array) {
	    //Undeclared array. Its size is inferred from the dimensions of
	    //the newly created node
	    vector<unsigned long> const &dim = node->dim();
	    for (unsigned int i = 0; i < dim.size(); ++i) {
		if (dim[i] == 0) {
		    CompileError(var, "Zero dimension for variable " +
				 var->name());
		}
	    }
	    symtab.addVariable(var->name(), dim);
	    array = symtab.getVariable(var->name());
	    array->insert(node, array->range());
	}
	else {
	    // Check if a node is already inserted into this range
	    SimpleRange range = VariableSubsetRange(var);
	    if (array->getSubset(range, _model)) {
		CompileError(var, "Attempt to redefine node",
			     var->name() + printRange(range));
	    }
	    array->insert(node, range);
	}
	_n_resolved++;
	_is_resolved[_n_relations] = true;
    }
    else if (_compiler_mode == CLEAN_UNRESOLVED) {
	/* 
	   Remove from the set of unresolved parameters, any array
	   subsets that are defined on the left hand side of a
	   relation
	*/
	ParseTree *var = rel->parameters()[0];
	SimpleRange range = VariableSubsetRange(var);
	_umap.erase(pair<string, Range>(var->name(), range));
	
	/*
	  In addition, the parameter might be a *subset* of a node
	  defined on the LHS of a relation.
	*/
	map<pair<string, Range>, set<unsigned long> >::iterator p = _umap.begin();
	while (p != _umap.end()) {
	    pair<string, Range> const &up = p->first;
	    if (up.first == var->name() && range.contains(up.second)) {
		_umap.erase(p++);
	    }
	    else {
		p++;
	    }
	}
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
    SimpleRange range = VariableSubsetRange(var);
    SimpleRange const &var_range = q->second.range();
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
       Called by traverseTree.

       For each node array, calculates the upper bounds of the subset
       expressions on the left-and side of all relations. The results
       are stored in the map _node_array_bounds which is then used to
       calculate the dimensions of the node arrays.
    */

    ParseTree const *var = p->parameters()[0];

    if (var->parameters().empty()) {
	//No index expession => No info on array size
	return;
    }

    SimpleRange new_range = VariableSubsetRange(var);
    vector<unsigned long> const &new_upper = new_range.last();

    string const &name = var->name();
    map<string, vector<unsigned long> >::iterator i = _node_array_bounds.find(name);
    if (i == _node_array_bounds.end()) {
	//Create a new entry
	_node_array_bounds[name] = new_upper;
    }
    else {
	//Check against the existing entry, and modify if necessary
	vector<unsigned long> & ubound = i->second;
	if (new_upper.size() != ubound.size()) {
	    CompileError(var, "Inconsistent dimensions for array", name);
	}
	else {
	    for (unsigned long j = 0; j < ubound.size(); ++j) {
		ubound[j] = max(ubound[j], new_upper[j]);
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

    _is_resolved = vector<bool>(_n_relations, false);
    for (unsigned long N = _n_relations; N > 0; N -= _n_resolved) {
	_n_resolved = 0;
	/* 
	   Here we use the abilitiy to sweep forwards and backwards
	   through the relations, allowing rapid compilation of models
	   written in both topological order and reverse topological
	   order. Without this facility, compilation of some large
	   models can be very slow.
	*/
	traverseTree(relations, &Compiler::allocate, true, true);
	if (_n_resolved == 0) {
	    break;
	}
    }
    _is_resolved.clear();
    
    if (_n_resolved == 0) {
	/*
	   Some nodes remain unresolved. We need to identify them and
	   print an informative error message for the user.
	   
	   The basic strategy is to identify nodes that appear on the
	   right-hand side of a relation but cannot be resolved.  Then
	   we eliminate the nodes that are defined on the left hand
	   side of a relation.  Whatever nodes are left are used in an
	   expression but never defined (either in a relation or in
	   the data).
	   
	   It may happen that all the unresolved nodes are defined on
	   the left hand side of a relation. This can happen if there
	   is a directed cycle (i.e. the nodes are defined in terms of
	   each other).
	*/

	
	/*
	  Step 2: Collect identifiers for unresolved parameters.
	  These will be held in the map _umap along with the line
	  numbers on which they are used.
	*/
	_compiler_mode = COLLECT_UNRESOLVED; //See getArraySubset
	traverseTree(relations, &Compiler::allocate);
	if (_umap.empty()) {
	    //Not clear what went wrong here, so throw a generic error message
	    throw runtime_error("Unable to resolve relations");
	}
	//Take a back-up copy of unresolved parameters
	map<pair<string,Range>, set<unsigned long> > umap_copy = _umap;

	/*
	  Step 3: Eliminate parameters that appear on the left of a relation
	*/
	_compiler_mode = CLEAN_UNRESOLVED;
	traverseTree(relations, &Compiler::allocate);

	/* 
	   Step 4: Informative error message for the user
	*/
	ostringstream oss;
	if (!_umap.empty()) {
	    //Undefined parameter message
	    oss << "Unable to resolve the following parameters:\n";
	    for (map<pair<string, Range>, set<unsigned long> >::iterator p =
		     _umap.begin(); p != _umap.end(); ++p)
	    {
		oss << p->first.first << printRange(p->first.second);
		oss << " (line ";
		set<unsigned long> const &lines = p->second;
		for (set<unsigned long>::const_iterator i = lines.begin();
		     i != lines.end(); ++i)
		{
		    if (i != lines.begin()) { oss << ", "; }
		    oss << *i;
		}
		oss << ")\n";
	    }
	    oss << "Either supply values for these nodes with the data\n"
		<< "or define them on the left hand side of a relation.";
	}
	else {
	    //Directed cycle message
	    oss << "Possible directed cycle involving some or all\n"
		<< "of the following nodes:\n";
	    for (map<pair<string, Range>, set<unsigned long> >::const_iterator p =
		     umap_copy.begin(); p != umap_copy.end(); ++p)
	    {
		oss << p->first.first << printRange(p->first.second) << "\n";
	    }
	}
	throw runtime_error(oss.str());
    }
}


void Compiler::traverseTree(ParseTree const *relations, CompilerMemFn fun,
			    bool resetcounter, bool reverse)
{
    /* 
       Traverse parse tree, expanding FOR loops and applying function
       fun to relations.

       In JAGS 4.3.0 we do a breadth first search, i.e. we process the
       relations within each block before expanding the for loops. This
       allows us to optionally reverse through the relations. See
       writeRelations for more details.
    */

    if (resetcounter) {
	_n_relations = 0;
    }

    vector<ParseTree*> const &relation_list = relations->parameters();
    for (vector<ParseTree*>::const_reverse_iterator p = relation_list.rbegin(); 
	 p != relation_list.rend(); ++p) 
    {
	TreeClass tc = (*p)->treeClass();
	if (tc == P_STOCHREL || tc == P_DETRMREL) {
	    (this->*fun)(*p);
	    _n_relations++;
	}
    }

    if (reverse) {

	unsigned int nrel = _n_relations; //Save current relation number
      
	for (vector<ParseTree*>::const_iterator p = relation_list.begin(); 
	     p != relation_list.end(); ++p) 
	{
	    //Reverse sweep through relations
	    TreeClass tc = (*p)->treeClass();
	    if (tc == P_STOCHREL || tc == P_DETRMREL) {
		_n_relations--;
		(this->*fun)(*p);
	    }
	}

	_n_relations = nrel; //Restore current relation number
    }
      
    for (vector<ParseTree*>::const_reverse_iterator p = relation_list.rbegin(); 
	 p != relation_list.rend(); ++p) 
    {
	if ((*p)->treeClass() == P_FOR) {
	    //Expand for loops
	    ParseTree *var = (*p)->parameters()[0];
	    if (!isNULL(CounterRange(var))) {
		Counter *counter = _countertab.pushCounter(var->name(),
							   CounterRange(var));
		for (; !counter->atEnd(); counter->next()) {
		    traverseTree((*p)->parameters()[1], fun, false, reverse);
		}
		_countertab.popCounter();
	    }
	}
    }
  
}


Compiler::Compiler(BUGSModel &model, map<string, SArray> const &data_table)
    : _model(model), _countertab(), 
      _data_table(data_table), _n_resolved(0), 
      _n_relations(0), _is_resolved(0), _compiler_mode(STANDARD),
      _index_expression(0), _index_nodes()
{
    if (_model.nodes().size() != 0)
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

  //checkDecNames(dec_list, _data_table);

  for (p = dec_list.begin() ; p != dec_list.end(); ++p) {
    ParseTree const *node_dec = *p;
    string const &name = node_dec->name();
    unsigned long ndim = node_dec->parameters().size();
    if (ndim == 0) {
	// Variable is scalar
	_model.symtab().addVariable(name, vector<unsigned long>(1,1));
    }
    else {
      // Variable is an array
	vector<unsigned long> dim(ndim);
	for (unsigned int i = 0; i < ndim; ++i) {
	    vector<unsigned long> dim_i;
	    if (!indexExpression(node_dec->parameters()[i], dim_i)) {
		CompileError(node_dec, "Unable to calculate dimensions of",
			     name);
	    }
	    if (dim_i.empty()) {
		CompileError(node_dec, "NULL dimension in declaration of", 
			     name);
	    }
	    if (dim_i.size() != 1) {
		CompileError(node_dec, 
			     "Vector-valued dimension in declaration of", name);
	    }
	    if (dim_i[0] <= 0) {
		CompileError(node_dec, "Non-positive dimension for node", name);
	    }
	    dim[i] = static_cast<unsigned int>(dim_i[0]);
	}
	_model.symtab().addVariable(name, dim);
    }

    //Check consistency with data, if supplied
    map<string, SArray>::const_iterator q = _data_table.find(name);
    if (q != _data_table.end()) {
	NodeArray const *array = _model.symtab().getVariable(name);
	if (q->second.range() != array->range()) {
	    string msg = string("Dimensions of ") + name + 
		" in declaration (" + printRange(array->range()) + 
		") conflict with dimensions in data (" + 
		printRange(q->second.range()) + ")";
	    CompileError(node_dec, msg);
	}
    }
  }
}

void Compiler::undeclaredVariables(ParseTree const *prelations)
{
    //Find parameters on RHS of relations that are never defined
    //checkVarNames(prelations, _data_table);
    
    // Get undeclared variables from data table
    for (map<string, SArray>::const_iterator p = _data_table.begin();
	 p != _data_table.end(); ++p) 
    {
	string const &name = p->first;
	NodeArray const *array = _model.symtab().getVariable(name);
	if (array) {
	    if (p->second.range() != array->range()) {
		//Should have been checked already in declareVariables
		throw logic_error("Dimension mismatch");
	    }
	}
	else {
	    _model.symtab().addVariable(name, p->second.dim(false));
	}
    }

    // Get names of variables on the LHS of a relation
    getLHSVars(prelations);

    // Infer the dimension of remaining nodes from the relations
    _compiler_mode = GET_DIMS;
    traverseTree(prelations, &Compiler::getArrayDim);
    _compiler_mode = STANDARD;
    
    map<string, vector<unsigned long> >::const_iterator i = _node_array_bounds.begin(); 
    for (; i != _node_array_bounds.end(); ++i) {
	if (_model.symtab().getVariable(i->first)) {
	    //Node already declared or defined by data. Check consistency 
	    NodeArray const * array = _model.symtab().getVariable(i->first);
	    vector<unsigned long> const &upper = array->range().upper();
	    if (upper.size() != i->second.size()) {
		string msg = "Dimension mismatch for variable ";
		msg.append(i->first);
		throw runtime_error(msg);
	    }
	    for (unsigned int j = 0; j < upper.size(); ++j) {
		if (i->second[j] <= 0 || i->second[j] > upper[j]) {
		    string msg =  string("Index out of range for variable ") + 
			i->first;
		    throw runtime_error(msg);
		}
	    } 
	}
	else {
	    //Variable not declared. Use inferred size
	    vector<unsigned long> const &upper = i->second;
	    unsigned long ndim = upper.size();
	    vector<unsigned long> dim(ndim);
	    for (unsigned int j = 0; j < ndim; ++j) {
		if (upper[j] <= 0) {
		    string msg = string("Invalid dimension for ") + i->first;
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

/* 
   We use construct-on-first-use for the lookup tables used by the
   compiler.  By dynamically allocating a table, we ensure that its
   destructor is never called - the memory is simply returned to the
   OS on exit.

   This fixes a nasty exit bug.  We cannot guarantee the order that
   static destructors are called in.  Therefore, a segfault can occur
   if a module tries to unload itself from a table that has already
   been destroyed.

   See also Model.cc, where the same technique is used for factory
   lists.
*/

DistTab &Compiler::distTab()
{
    static DistTab *_disttab = new DistTab();
    return *_disttab;
}

FuncTab &Compiler::funcTab()
{
    static FuncTab *_functab = new FuncTab();
    return *_functab;
}

ObsFuncTab &Compiler::obsFuncTab()
{
    static ObsFuncTab *_oftab = new ObsFuncTab();
    return *_oftab;
}

MixtureFactory& Compiler::mixtureFactory1()
{
   return _mixfactory1;
}

MixtureFactory& Compiler::mixtureFactory2()
{
   return _mixfactory2;
}

BUGSModel &Compiler::model() const
{
    return _model;
}

} //namespace jags
