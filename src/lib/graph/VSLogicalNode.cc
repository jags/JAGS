#include <config.h>
#include <graph/VSLogicalNode.h>
#include <function/ScalarFunction.h>

#include <stdexcept>
#include <vector>
#include <math.h>

using std::vector;
using std::logic_error;

namespace jags {

static vector<unsigned long> mkDim(vector<Node const *> const &parameters)
{
    vector<unsigned long> dim(1,1);
    bool scalar = true;
    for (unsigned long i = 0; i < parameters.size(); ++i) {
	if (parameters[i]->length() > 1) {
	    if (scalar) {
		dim = parameters[i]->dim();
		scalar = false;
	    }
	    else if (dim != parameters[i]->dim()) {
		throw logic_error("Incompatible dimensions in VSLogicalNode");
	    }
	}
    }
    return dim;
}

static vector<bool> mkIsVector(vector<Node const *> const &parameters)
{
    vector<bool> ans(parameters.size());
    for (unsigned long i = 0; i < parameters.size(); ++i) {
	ans[i] = (parameters[i]->length() > 1);
    }
    return ans;
}

VSLogicalNode::VSLogicalNode(ScalarFunction const *function,
			     unsigned int nchain,
			     vector<Node const *> const &parameters)
    : LogicalNode(mkDim(parameters), nchain, parameters, function),
      _func(function), _isvector(mkIsVector(parameters))
{
    if (isFixed()) {
	for (unsigned int ch = 0; ch < nchain; ++ch) {
	    deterministicSample(ch);
	}
    }
}

void VSLogicalNode::deterministicSample(unsigned int chain)
{
    double *ans = _data + chain * _length;
    vector<double const *> par(_parameters[chain]);
	
    for (unsigned int i = 0; i < _length; ++i) {
	ans[i] = _func->evaluate(par);
	for (unsigned int j = 0; j < par.size(); ++j) {
	    if (_isvector[j])
		++par[j];
	}
    }
}

bool VSLogicalNode::checkParentValues(unsigned int chain) const
{
    vector<double const *> par(_parameters[chain]);

    for (unsigned int i = 0; i < _length; ++i) {
	if (!_func->checkParameterValue(par))
	    return false;
	for (unsigned long j = 0; j < par.size(); ++j) {
	    if (_isvector[j])
		++par[j];
	}
    }
    return true;
}

    void VSLogicalNode::gradient(double *grad, Node const *arg,
				 unsigned int chain) const
    {
	vector<double const *> param(_parameters[chain]);

	auto par = parents();
	for (unsigned long i = 0; i < par.size(); ++i) {
	    if (par[i] != arg) continue;
	    
	    for (unsigned int l = 0; l < _length; ++l) {
		grad[l] += _func->gradient(param, i);
		for (unsigned int k = 0; k < param.size(); ++k) {
		    if (_isvector[k]) ++param[k];
		}
	    }  
	}
    }

    
    /*
DeterministicNode *
VSLogicalNode::clone(vector<Node const*> const &parents) const
{
    return new VSLogicalNode(_func, parents);
}
    */

} //namespace jags
