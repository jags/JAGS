#include <config.h>
#include <graph/ArrayLogicalNode.h>
#include <function/FuncError.h>
#include <function/ArrayFunction.h>
#include <util/dim.h>

#include <stdexcept>
#include <vector>
#include <string>

using std::vector;
using std::string;
using std::logic_error;

namespace jags {

static vector<unsigned int> mkDim(ArrayFunction const *func, 
				  vector<Node const *> const &parents)
{
    /* 
       Calculates dimension of logical node as a function of its
       parameters.
    */

    unsigned long N = parents.size();
    vector<vector<unsigned int> > dims(N);
    vector<bool> fixed(N), discrete(N);
    vector<double const *> values(N);
    for (unsigned int j = 0; j < N; ++j) {
	dims[j] = parents[j]->dim();
	fixed[j] = parents[j]->isFixed();
	values[j] = parents[j]->value(0);
	discrete[j] = parents[j]->isDiscreteValued();
    }
    if (!func) {
	throw logic_error("NULL function in ArrayLogicalNode constructor");
    }
    if (!checkNPar(func, dims.size())) {
	throw FuncError(func, "Incorrect number of parameters");
    }
    if (!func->checkParameterDim(dims)) {
	throw FuncError(func, "Non-conforming parameters");
    }
    if (!func->checkParameterFixed(fixed)) {
	throw FuncError(func, "Expected parameters with fixed values");
    }
    if (!func->checkParameterDiscrete(discrete)) {
	throw FuncError(func, "Failed check for discrete-valued parameters");
    }
    return func->dim(dims, values);
}

static vector<vector<unsigned int> > const &
mkParameterDims(vector<Node const *> const &parameters) {
    vector<vector<unsigned int> > dims(parameters.size());
    for (unsigned int j = 0; j < parameters.size(); ++j) {
        dims[j] = parameters[j]->dim();
    }
    return getUnique(dims);
}

ArrayLogicalNode::ArrayLogicalNode(ArrayFunction const *function,
				   unsigned int nchain,
				   vector<Node const *> const &parameters)
    : LogicalNode(mkDim(function,parameters), nchain, parameters, function),
      _func(function), _dims(mkParameterDims(parameters))
{
    initializeFixed();
}

void ArrayLogicalNode::deterministicSample(unsigned int chain)
{
    _func->evaluate(_data + chain * _length, _parameters[chain], _dims);
}

bool ArrayLogicalNode::checkParentValues(unsigned int chain) const
{
    return _func->checkParameterValue(_parameters[chain], _dims);
}

    /*
DeterministicNode *
ArrayLogicalNode::clone(vector<Node const*> const &parents) const
{
    return new ArrayLogicalNode(_func, parents);
}
    */
    
} //namespace jags
