#include <config.h>

#include <compiler/Compiler.h>
#include <compiler/ParseTree.h>
#include <model/SymTab.h>
#include <graph/LogicalNode.h>
#include <graph/ConstantNode.h>
#include <graph/StochasticNode.h>
#include <graph/DevianceNode.h>
#include <graph/Graph.h>
#include <sarray/RangeIterator.h>

#include "MixCompiler.h"

#include <utility>
#include <vector>
#include <stdexcept>
#include <cfloat>
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


/* FIXME: Use this everywhere */
static long asInteger(double fval)
{
  if (fval >= LONG_MAX || fval <= LONG_MIN) {
    throw runtime_error("double value out of range for conversion to long");
  }
  long ival;
  if (fval > 0) {
    ival = (long) (fval + DBL_EPSILON);
  }
  else {
    ival = (long) (fval - DBL_EPSILON);
  }
  if (fabs(fval - ival) > DBL_EPSILON) {
    throw runtime_error("Invalid integer conversion");
  }
  return ival;
}

double Compiler::constFromTable(ParseTree const *p)
{
  // Try evaluating constant expression from data table
  map<string,SArray>::const_iterator i = _data_table.find(p->name());
  if (i == _data_table.end()) {
    return JAGS_NA;
  }
  Range range = getRange(p->parameters(), i->second.range());
  if (isNULL(range)) {
    return JAGS_NA;
  }
  else {
    // Range expression successfully evaluated
    if (range.length() > 1) {
      throw runtime_error("Vector value in constant expression");
    }
    long offset = i->second.range().leftOffset(range.lower());
    return i->second.value()[offset];
  }
}

double Compiler::constFromNode(ParseTree const*p)
{
  //Evaluate constant expression if it corresponds to a node with
  //a fixed value

  NodeArray *array = _symtab.getVariable(p->name());
  if (array) {
    // We can't call getConstantRange() because we aren't sure that
    // the range expression can be evaluated.
    Range subset_range = getRange(p->parameters(), array->range());
    if (isNULL(subset_range)) {
      return JAGS_NA;
    }
    else if (subset_range.length() > 1) {
      throw runtime_error(string("Vector value ") + p->name() +
			  print(subset_range) + " in constant expression");
    }
    else {
      Node *node = array->getSubset(subset_range);
      if (node && node->isObserved()) {
	return *node->value(0);
      }
      else {
	return JAGS_NA;
      }
    }
  }
  else {
    return JAGS_NA;
  }
}

bool Compiler::constantExpression(ParseTree const *p, double &value)
{
  /* 
     Try to evaluate a constant expression.  An expression can
     only be evaluated if it is a function of
     1) Counters, or
     2) a value given in the data table, or
     3) a fixed Node.
     Functions are not allowed, but inline operators ("*", "/", "+", "-") 
     are.

     If the expression can be evaluated, then we return true and the
     expression value is written to the value argument. If it cannot
     be evaluated (because it depends on missing data) then we return
     false.
  */
  Counter *counter;
  double arg1, arg2;
  std::vector<ParseTree*> const &parameters = p->parameters();

  switch (p->treeClass()) {
  case P_VAR: 
    // Is it a counter?
    counter = _countertab.getCounter(p->name());
    if (counter) {
      value = (*counter)[0];
      return true;
    }
    else {
      // Is it a variable? Try reading value from the data table
      value = constFromTable(p);
      if (value != JAGS_NA) {
	return true;
      }
      else {
	// Failing that, see if it is a node with a fixed value
	value = constFromNode(p);
	return (value != JAGS_NA);
      }
    }
    break;
  case P_VALUE:
    value =  p->value();
    return true;
    break;
  case P_OPERATOR:
    switch(p->getOperator()) {
    case OP_ADD:
      if (parameters.size() < 2) {
	throw logic_error("ADD must have at least two arguments in constant expression");
      }
      arg1 = 0;
      for (unsigned int i = 0; i < parameters.size(); ++i) {
	if (!constantExpression(parameters[i], arg2)) {
	  return false;
	}
	arg1 += arg2;
      }
      value = arg1;
      return true;
      break;
    case OP_SUBTRACT:
      if (parameters.size() != 2) {
	throw logic_error("SUBTRACT must have two arguments in constant expression");
      }
      if (!constantExpression(parameters[0], arg1)) {
	return false;
      }
      if (!constantExpression(parameters[1], arg2)) {
	return false;
      }
      value = arg1 - arg2;
      return true;
      break;
    case OP_MULTIPLY:
      if (parameters.size() < 2) {
	throw logic_error("MULTIPLY must have at least two arguments in constant expression");
      }
      arg1 = 1;
      for (unsigned int i = 0; i < parameters.size(); ++i) {
	if (!constantExpression(parameters[i], arg2)) {
	  return false;
	}
	arg1 *= arg2;
      }
      value = arg1;
      return true;
      break;
    case OP_DIVIDE:
      if (parameters.size() != 2) {
	throw logic_error("DIVIDE must have two arguments in constant expression");
      }
      if (!constantExpression(parameters[0], arg1)) {
	return false;
      }
      if (!constantExpression(parameters[1], arg2)) {
	return false;
      }
      value = arg1/arg2;
      return true;
      break;
    case OP_NEG:
      if (parameters.size() != 1) {
	throw logic_error("NEG must have 1 argument in constant expression");
      }
      if (!constantExpression(parameters[0], arg1)) {
	return false;
      }
      value = -arg1;
      return true;
      break;
    case OP_NONE:
      throw logic_error("Bad constant expression");
      break;
    }
    return false;
    break;
  default:
    throw logic_error("Expected variable, value or expression");
  }
}

bool Compiler::indexExpression(ParseTree const *p, long &value)
{
  /* 
     Evaluate an index expression.
     
     Index expressions occur in three contexts:
     1) In the limits of a "for" loop
     2) On the left hand side of a relation
     3) On the right hand side of a relation
     
     In cases 1) and 2) we expect to be able to evaluate the index
     expression, and an exception will be thrown if it cannot be
     evaluated (strict == true). In case 3) it is not an error if an
     index expression cannot be evaluated (strict == false), but we try
     anyway for efficiency reasons.
  */

  double fvalue;
  if (!constantExpression(p, fvalue)) {
    return false;
  }
  else {
    value = asInteger(fvalue);
    return true;
  }
}

Range Compiler::getRange(vector<ParseTree*> const &range_list, 
			 Range const &default_range)
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
  
  if (range_list.empty()) {
    /* An empty range expression implies the default range, if it exists,
       or a scalar value, if it does not */
    if (isNULL(default_range))
      return Index(1);
    else 
      return default_range;
  }

  // Check size and integrity of range expression
  unsigned int size = range_list.size();
  if (!isNULL(default_range) && size != default_range.ndim(false)) {
    throw logic_error("Default range does not match dimension of range expression");
  }
  for (unsigned int i = 0; i < size; ++i) {
    if (range_list[i]->treeClass() != P_RANGE) {
      throw logic_error("Malformed parse tree. Expected range expression");
    }
  }
  
  // Now step through and evaluate lower and upper index expressions
  Index lower(size), upper(size);
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
  NodeArray *array = _symtab.getVariable(name);
  if (array) {
    // It's a declared node
    vector<ParseTree*> const &range_list = var->parameters();
    
    bool ok = true;
    if (range_list.empty()) {
      ok = (array->range().ndim(false) == 1);
    }
    else {
      ok = (range_list.size() == array->range().ndim(false));
    }
    if (!ok) {
      throw runtime_error(string("Dimension mismatch in subset expression ")
			  + "of variable " + name);
    }
    Range range = getRange(var->parameters(), array->range());
    if (isNULL(range)) {
      throw runtime_error(string("Missing values in subset expression ") 
			  + "of variable " + name);
    }
    return range;
  }
  else {
    // Undeclared node
    Range range = getRange(var->parameters(), Range());
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
  long lower;
  if(!indexExpression(prange->parameters()[0], lower)) {
    throw runtime_error(string("Unable to evaluate lower index of counter ")
			+ var->name());
  }
  long upper;
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
    Index ind_lower(1), ind_upper(1);
    ind_lower[0] = lower;
    ind_upper[0] = upper;
    
    return Range(ind_lower, ind_upper);
  }
}

Node* Compiler::VarGetNode(ParseTree const *var)
{
  if (var->treeClass() != P_VAR) {
    throw logic_error("Expecting variable expression");
  }
  NodeArray *array = _symtab.getVariable(var->name());
  Range range = VariableSubsetRange(var);
  Node *node = array->find(range);
  
  return node;
}

Node * Compiler::getSubSetNode(ParseTree const *var)
{
  if (var->treeClass() != P_VAR) {
    throw logic_error("Expecting variable expression");
  }

  NodeArray *array = _symtab.getVariable(var->name());
  if (array == 0) {
    throw runtime_error(string("Unknown variable ") + var->name());
  }
  Range subset_range = getRange(var->parameters(), array->range());
  if (isNULL(subset_range)) {
    return 0;
  }
  if (!array->range().contains(subset_range)) {
    throw runtime_error(string("Subset ") + var->name() + print(subset_range)
			+ " out of range");
  }
  Node *node = array->getSubset(subset_range);
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

  switch(p->treeClass()) {
  case P_VALUE:
    node = _constantfactory.getConstantNode(p->value());
    break;
  case P_VAR:
    {
      Counter *counter = _countertab.getCounter(p->name()); //A counter
      if (counter) {
	node = _constantfactory.getConstantNode((*counter)[0]);
      }
      else {
	NodeArray *array = _symtab.getVariable(p->name());
	if (array == 0) {
	  throw runtime_error(string("Unknown parameter ") + p->name());
	}
	Range subset_range = getRange(p->parameters(), array->range());
	if (isNULL(subset_range)) {
	  node = getMixtureNode(p, this); //A stochastic subset
	} 
	else {
	  node = getSubSetNode(p); //A fixed subset
	} 
      }
    }
    break;
  default:
    throw logic_error("Expecting value or variable expression");
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
  Function const *func = 0;
  
  switch (t->treeClass()) {
  case P_FUNCTION:
    func = functab.find(t->name());
    if (func == 0) {
      string msg("Unable to find function ");
      msg.append(t->name());
      throw runtime_error(msg);
    }    
    break;
  case P_OPERATOR:
    switch(t->getOperator()) {
    case OP_ADD:
      func = functab.find("+");
      break;
    case OP_SUBTRACT:
      func = functab.find("-");
      break;
    case OP_MULTIPLY:
      func = functab.find("*");
      break;
    case OP_DIVIDE:
      func = functab.find("/");
      break;
    case OP_NEG:
      func = functab.find("NEG");
      break;
    case OP_NONE:
      throw logic_error("Bad operator expression");
      break;
    }
    if (func == 0) {
      throw logic_error("Unable to find operator");
    }
    break;
  default:
    throw logic_error("Malformed parse tree: Expected expression");
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

//debuggin
#include <iostream>
Node* Compiler::getParameter(ParseTree const *t)
{
    vector<Node*> parents;
    Node *node = 0;

    switch (t->treeClass()) {
    case P_VALUE:
	node =  _constantfactory.getConstantNode(t->value());
	break;
    case P_VAR:
	node = getArraySubset(t);
	break;
    case P_FUNCTION: case P_OPERATOR:
	if (getLogicalParameterVector(t, parents)) {
	    node = _logicalfactory.getLogicalNode(getFunction(t, funcTab()), parents);
	}
	break;
    default:
	throw  logic_error("Malformed parse tree. Expected value, variable or expression");
	break;
    }

    /* Initialize deterministic nodes now, if they are functions of 
       data, to aid in node recycling */
    if (node && !node->isVariable()) {
      node->initialize();
    }
    return node;
}

bool Compiler::getLogicalParameterVector(ParseTree const *t,
					 vector<Node*> &parents)
{
  if (!parents.empty()) {
    throw logic_error("parent vector must be empty in getLogicalParameterVector");
  }

  switch (t->treeClass()) {
  case P_FUNCTION: case P_LINK: case P_OPERATOR:
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
    throw logic_error("Invalid Parse Tree.  Expected function or operator.");
  }
  return true;
}

Node * Compiler::allocateStochastic(ParseTree const *stoch_relation)
{
  ParseTree const *distribution = stoch_relation->parameters()[1];  
  vector<Node*> parameters;
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
      lBound = getArraySubset(ll);
      if (!lBound) {
	return 0;
      }
    }
    if (ul) {
      uBound = getArraySubset(ul);
      if (!uBound) {
	return 0;
      }
    }
  }

  return new StochasticNode(getDistribution(stoch_relation, distTab()),
			    parameters, lBound, uBound);
}

Node * Compiler::allocateLogical(ParseTree const *rel)
{
  ParseTree *expression = rel->parameters()[1];
  Node *node = 0;
  vector <Node*> parents;

  switch (expression->treeClass()) {
  case P_VALUE: 
    node = new ConstantNode(expression->value(), _symtab.nchain());
    _graph.add(node); //FIXME: Why are we adding here?
    break;
  case P_VAR: case P_FUNCTION: case P_OPERATOR:
    node = getParameter(expression);
    break;
  case P_LINK:
    if (getLogicalParameterVector(expression, parents)) {
      node = _logicalfactory.getLogicalNode(getLink(expression, funcTab()), 
					    parents);
    }
    break;
  default:
    throw logic_error("Malformed parse tree in Compiler::allocateLogical");
  }

  return node;
}

#include <iostream>
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

  if (node) {
    /* Check if a node is already inserted into this range */
    ParseTree *var = rel->parameters()[0];
    NodeArray *array = _symtab.getVariable(var->name());
    Range range = VariableSubsetRange(var);
    if (array->find(range)) {
      throw runtime_error(string("Attempt to redefine node ") +
			  var->name() + print(range));
    }
    array->insert(node, range);
    _n_resolved++;
    _is_resolved[_n_relations] = true;
  }
}

#include <iostream>
void Compiler::setConstantMask(ParseTree const *rel)
{
  ParseTree const *var = rel->parameters()[0];
  string const &name = var->name();
  Range range = VariableSubsetRange(var);
  map<string,vector<bool> >::iterator p = _constant_mask.find(name);
  if (p == _constant_mask.end()) {
    return;
  }
  map<string,SArray>::const_iterator q = _data_table.find(name);
  if (q == _data_table.end()) {
    throw logic_error ("Error in Compiler::setConstantMask");
  }
  Range const &var_range = q->second.range();
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

  Range new_range = VariableSubsetRange(var);

  map<string, vector<Index> >::iterator i = _node_array_ranges.find(name);
  if (i == _node_array_ranges.end()) {
    //Create a new entry
    vector<Index> ivec;
    ivec.push_back(new_range.lower());
    ivec.push_back(new_range.upper());
    _node_array_ranges.insert(pair<const string, vector<Index> >(name,ivec));
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
    pair<string,vector<bool> > apair;
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

  _symtab.writeData(temp_data_table);
}

void Compiler::writeVariableData()
{
  //Create a temporary copy of the data table with no data for
  //constant nodes
  map<string, SArray> temp_data_table = _data_table;
  map<string, SArray>::iterator p;
  for(p = temp_data_table.begin(); p != temp_data_table.end(); ++p) {
    string const &name = p->first;
    SArray &temp_data = p->second;
    vector<bool> const &mask = _constant_mask.find(name)->second;
    for (unsigned long i = 0; i < temp_data.length(); ++i) {
      if (mask[i]) {
	temp_data.setValue(JAGS_NA, i);
      }
    }
  }

  _symtab.writeData(temp_data_table);
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

  writeVariableData();
  collectNodes();
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

Compiler::Compiler(Graph &graph, SymTab &symtab,
		   map<string, SArray> const &data_table)
    : _graph(graph), _symtab(symtab), _countertab(), 
      _data_table(data_table), _n_resolved(0), 
      _n_relations(0), _is_resolved(0), _strict_resolution(false),
      _constantfactory(symtab.nchain())
{
  if (graph.size() != 0)
    throw invalid_argument("Non empty graph in Compiler constructor");
  if (symtab.size() != 0)
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
      _symtab.addVariable(name);
    }
    else {
      // Variable is an array
      Index dim(ndim);
      for (unsigned int i = 0; i < ndim; ++i) {
	if (!indexExpression(node_dec->parameters()[i], dim[i])) {
	  throw runtime_error(string("Unable to calculate dimensions of node ")
			      + name);
	}
      }
      _symtab.addVariable(name, dim);
    }
  }
}

void Compiler::collectNodes()
{
  vector<Node*> nodes;
  _symtab.getNodes(nodes);
  _logicalfactory.graph().getNodes(nodes);
  _constantfactory.graph().getNodes(nodes);
  _mixfactory.graph().getNodes(nodes);

  for (vector<Node*>::iterator p = nodes.begin(); p != nodes.end(); ++p) {
    _graph.add(*p);
  }
}

void Compiler::undeclaredVariables(ParseTree const *prelations)
{
  // Get undeclared variables from data table
  map<string, SArray>::const_iterator p = _data_table.begin();
  for (; p != _data_table.end(); ++p) {
    string const &name = p->first;
    if (!_symtab.getVariable(name)) {
      _symtab.addVariable(name, p->second.dim(false));
    }
  }

  // Infer the dimension of remaining nodes from the relations
  traverseTree(prelations, &Compiler::getArrayDim);
  map<string, vector<Index> >::const_iterator i = _node_array_ranges.begin(); 
  for (; i != _node_array_ranges.end(); ++i) {
    if (_symtab.getVariable(i->first)) {
       //Node already declared. Check consistency 
       NodeArray const * array = _symtab.getVariable(i->first);
       Index const &upper = array->range().upper();
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
       _symtab.addVariable(i->first, i->second[1]);
    }
  }
  //_node_array_ranges.clear(); (Need to keep this now)
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

SymTab &Compiler::symTab() const
{
   return _symtab;
}

Graph &Compiler::graph() const
{
  return _graph;
}
