#include <config.h>
#include <graph/LogicalNode.h>
#include <function/Function.h>
#include <graph/NodeError.h>
#include <graph/NodeNameTab.h>
#include <sarray/Index.h>

#include <stdexcept>
#include <vector>
#include <string>
#include <math.h>

using std::vector;
using std::string;
using std::set;
using std::logic_error;
using std::runtime_error;

#include <iostream>
static Index mkDim(Function const *func, vector<Node*> const &parents)
{
  /* 
     Calculates dimension of logical node as a function of its
     parameters.
  */

  vector<Index> parameter_dims(parents.size());
  for (unsigned long j = 0; j < parents.size(); ++j) {
    parameter_dims[j] = parents[j]->dim(true);
  }
  if (!func->checkParameterLength(parameter_dims.size())) {
    //FIXME: logic_error or runtime_error?
    throw runtime_error(string("Incorrect number of parameters for function ")
		      + func->name());
  }
  if (!func->checkParameterDim(parameter_dims)) {
    throw runtime_error(string("Non-conforming parameters for function ")
		      + func->name());
  }
  return func->dim(parameter_dims);
}
  
LogicalNode::LogicalNode(Function const *function, 
			 vector<Node*> const &parameters)
  : DeterministicNode(mkDim(function,parameters), parameters),
    _func(function),
    _parameters(nchain())
{
  
  for (unsigned int n = 0; n < nchain(); ++n) {
    _parameters[n].reserve(parameters.size());
    for (unsigned long j = 0; j < parameters.size(); ++j) {
      _parameters[n].push_back(parameters[j]->data(n));
    }
  }

  vector<bool> mask(parameters.size());
  for (unsigned long j = 0; j < parameters.size(); ++j) {
    mask[j] = parameters[j]->isDiscreteValued();
  }
  if (_func->isDiscreteValued(mask)) {
     setDiscreteValued();
  }
}

LogicalNode::~LogicalNode()
{
}

string LogicalNode::name(NodeNameTab const &name_table) const
{
  string name = name_table.getName(this);
  if (!name.empty())
    return name;

  vector<Node *> const &par = parents();
  string const &fname = _func->name();
  if (fname == "*" || fname == "/" || fname == "+" || fname == "-")
    {
      name.append("(");
      for (unsigned int i = 0; i < par.size(); ++i) {
	if (i > 0) {
	  name.append(fname);
	}
	name.append(par[i]->name(name_table));
      }
      name.append(")");
    }
  else if (fname == "NEG") 
    {
      name.append("(-");
      name.append(par[0]->name(name_table));
      name.append(")");
    }
  else
    {
      name.append(fname);
      name.append("(");
      for (unsigned int i = 0; i < par.size(); ++i) {
	if (i > 0) {
	  name.append(",");
	}
	name.append(par[i]->name(name_table));
      }
      name.append(")");
    }

  return name;
}

vector<SArray const *> const &LogicalNode::parameters(unsigned int chain) const
{
  return _parameters[chain];
}

Function const* LogicalNode::function() const
{
  return _func;
}

void LogicalNode::deterministicSample(unsigned int chain)
{
  if (!_func->checkParameterValue(_parameters[chain])) {
    throw NodeError(this, "Invalid parameter values for LogicalNode");
  }
  _func->evaluate(*_data[chain], _parameters[chain]);
}

LogicalNode const *asLogical(Node const *node)
{
  return dynamic_cast<LogicalNode const*>(node);
}

bool isLogical(Node const *node)
{
  return dynamic_cast<LogicalNode const*>(node);
}

bool LogicalNode::isLinear(std::set<Node*> const &parameters, bool fixed) const
{
    vector<bool> mask(parents().size());
    for (unsigned int i = 0; i < parents().size(); ++i) {
	mask[i] = parameters.count(parents()[i]);
    }

    vector<bool> fixed_mask;
    if (fixed) {
	for (unsigned int i = 0; i < parents().size(); ++i) {
	    fixed_mask.push_back(parents()[i]->isObserved());
	}
    }

    return _func->isLinear(mask, fixed_mask);
}

bool LogicalNode::isScale(std::set<Node*> const &parameters, bool fixed) const
{
    unsigned int nparam = 0;
    unsigned int index = 0;
    for (unsigned int i = 0; i < parents().size(); ++i) {
	if (parameters.count(parents()[i])) {
	    nparam++;
	    index = i;
	}
    }
    if (nparam == 0)
	return true;
    else if (nparam > 1)
	return false;
    
    vector<bool> fixed_mask;
    if (fixed) {
	for (unsigned int i = 0; i < parents().size(); ++i) {
	    fixed_mask.push_back(parents()[i]->isObserved());
	}
    }

    return _func->isScale(index, fixed_mask);
}
