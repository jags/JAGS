#include <config.h>
#include <graph/ScalarLogicalNode.h>
#include <graph/NodeError.h>
#include <function/ScalarFunction.h>
#include <graph/GraphMarks.h>
#include <graph/Graph.h>
#include <util/dim.h>

#include <stdexcept>
#include <vector>
#include <string>
#include <math.h>

using std::vector;
using std::string;
using std::set;
using std::logic_error;

namespace jags {

ScalarLogicalNode::ScalarLogicalNode(ScalarFunction const *function,
				     unsigned int nchain,
				     vector<Node const *> const &parameters)
    : LogicalNode(vector<unsigned long>(1,1), nchain, parameters, function),
      _func(function)
{
    if (!function) {
	throw logic_error("NULL function in ScalarLogicalNode constructor");
    }
    for (unsigned long j = 0; j < parameters.size(); ++j) {
	if (isFlat(parameters[j]->dim())) {
	    string msg("Invalid zero-length parameter to function ");
	    msg.append(function->name());
	    throw NodeError(parameters[j], msg);
	}
	else if (!isScalar(parameters[j]->dim())) {
	    string msg("Invalid non-scalar parameter to function ");
	    msg.append(function->name());
	    throw NodeError(parameters[j], msg);
	}
    }

    initializeFixed();
}

void ScalarLogicalNode::deterministicSample(unsigned int chain)
{
    _data[chain] = _func->evaluate(_parameters[chain]);
}

bool ScalarLogicalNode::checkParentValues(unsigned int chain) const
{
    return _func->checkParameterValue(_parameters[chain]);
}

    /*
DeterministicNode *
ScalarLogicalNode::clone(vector<Node const*> const &parents) const
{
    return new ScalarLogicalNode(_func, parents);
}
    */
    
} //namespace jags
