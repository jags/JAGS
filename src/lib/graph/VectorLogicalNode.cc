#include <config.h>
#include <graph/VectorLogicalNode.h>
#include <function/FuncError.h>
#include <function/VectorFunction.h>
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

static unsigned long valueLength(VectorFunction const *func, 
				vector<Node const *> const &parents)
{
    /* 
       Calculates length of vector logical node as a function of its
       parents.
    */
    
    unsigned long N = parents.size();
    vector<unsigned long> lengths(N);
    vector<bool> fixed(N), discrete(N);
    vector<double const *> values(N);
    for (unsigned long j = 0; j < N; ++j) {
	lengths[j] = parents[j]->length();
	fixed[j] = parents[j]->isFixed();
	values[j] = parents[j]->value(0);
	discrete[j] = parents[j]->isDiscreteValued();
    }
    
    if (!func) {
	throw logic_error("NULL function in VectorLogicalNode constructor");
    }
    if (!checkNPar(func, N)) {
	throw FuncError(func, "Incorrect number of parameters");
    }
    if (!func->checkParameterLength(lengths)) {
	throw FuncError(func, "Non-conforming parameters");
    }
    if (!func->checkParameterFixed(fixed)) {
	throw FuncError(func, "Expected parameters with fixed values");
    }
    if (!func->checkParameterDiscrete(discrete)) {
	throw FuncError(func, "Failed check for discrete-valued parameters");
    }
    return func->length(lengths, values);
}

static vector<unsigned long> const &
parameterLengths(vector<Node const *> const &parameters) {
    vector<unsigned long> lengths(parameters.size());
    for (unsigned long j = 0; j < parameters.size(); ++j) {
        lengths[j] = parameters[j]->length();
	
    }
    return getUnique(lengths);
}

VectorLogicalNode::VectorLogicalNode(VectorFunction const *function,
				     unsigned int nchain,
				     vector<Node const *> const &parameters)
    : LogicalNode(vector<unsigned long>(1,valueLength(function, parameters)),
		  nchain, parameters, function),
      _func(function), _lengths(parameterLengths(parameters))
{
    initializeFixed();
}

void VectorLogicalNode::deterministicSample(unsigned int chain)
{
    _func->evaluate(_data + chain * _length, _parameters[chain], _lengths);
}

bool VectorLogicalNode::checkParentValues(unsigned int chain) const
{
    return _func->checkParameterValue(_parameters[chain], _lengths);
}

    void VectorLogicalNode::gradient(double *grad, Node const *arg,
				     unsigned int chain) const
    {
	auto par = parents();
	for (unsigned int i = 0; i < par.size(); ++i) {
	    if (par[i] == arg) {
		_func->gradient(grad, _parameters[chain], _lengths, i);
	    }
	}
    }

    /*
DeterministicNode *
VectorLogicalNode::clone(vector<Node const*> const &parents) const
{
    return new VectorLogicalNode(_func, parents);
}
    */
    
} //namespace jags
