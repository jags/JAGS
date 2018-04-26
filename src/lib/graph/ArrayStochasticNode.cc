#include <config.h>
#include <graph/ArrayStochasticNode.h>
#include <distribution/DistError.h>
#include <distribution/ArrayDist.h>
#include <util/nainf.h>
#include <util/dim.h>

#include <vector>
#include <string>
#include <algorithm>

using std::vector;
using std::string;
using std::copy;
using std::min;
using std::max;

namespace jags {

static vector<unsigned long> mkDim(ArrayDist const *dist, 
				  vector<Node const *> const &parents)
{
    /* 
       Calculates dimension of array stochastic node as a function of its
       parents
    */

    if (!checkNPar(dist, parents.size())) {
	throw DistError(dist, "Incorrect number of parameters");
    }
    vector<vector<unsigned long> > parameter_dims(parents.size());
    for (unsigned long j = 0; j < parents.size(); ++j) {
	parameter_dims[j] = parents[j]->dim();
    }
    if (!dist->checkParameterDim(parameter_dims)) {
	throw DistError(dist, "Non-conforming parameters");
    }
    return dist->dim(parameter_dims);
}

static vector<vector<unsigned long> > const &
mkParameterDims(vector<Node const *> const &parameters) {
    vector<vector<unsigned long> > dims(parameters.size());
    for (unsigned long j = 0; j < parameters.size(); ++j) {
        dims[j] = parameters[j]->dim();
    }
    return getUnique(dims);
}

ArrayStochasticNode::ArrayStochasticNode(ArrayDist const *dist,
					 unsigned int nchain,
					 vector<Node const *> const &params,
					 Node const *lower, Node const *upper,
					 double const *data,
					 unsigned long length)
    : StochasticNode(mkDim(dist, params), nchain, dist, params, lower, upper),
      _dist(dist), _dims(mkParameterDims(params))
{
    if (!dist->checkParameterDim(_dims)) {
	throw DistError(dist, "Invalid parameter dimensions");
    }
}

double ArrayStochasticNode::logDensity(unsigned int chain, PDFType type) const
{
    if(!_dist->checkParameterValue(_parameters[chain], _dims))
	return JAGS_NEGINF;
    
    return _dist->logDensity(_data + _length * chain, _length, type,
			     _parameters[chain], _dims,
			     lowerLimit(chain), upperLimit(chain));
}

void ArrayStochasticNode::randomSample(RNG *rng, unsigned int chain)
{
    _dist->randomSample(_data + _length * chain, _length,
			_parameters[chain], _dims, 
			lowerLimit(chain), upperLimit(chain), rng);
}  

void ArrayStochasticNode::truncatedSample(RNG *rng, unsigned int chain,
					  double const *lower, 
					  double const *upper)
{
    /*
      Some complexity required to deal with case where node is already
      truncated. Note that this code is identical to that in
      VectorStochasticNode
    */
    
    double const *l = lowerLimit(chain);
    double *lv = 0;

    if (l || lower) {
	lv = new double[_length];
	if (l && lower) {
	    for (unsigned long i = 0; i < _length; ++i) {
		lv[i] = min(l[i], lower[i]);
	    }
	}
	else if (l) {
	    copy(l, l + _length, lv);
	}
	else if (lower) {
	    copy(lower, lower + _length, lv);
	}
    }

    double const *u = upperLimit(chain);
    double *uv = 0;
    if (u || upper) {
	uv = new double[_length];
	if (u && upper) {
	    for (unsigned long i = 0; i < _length; ++i) {
		uv[i] = max(u[i], upper[i]);
	    }
	}
	else if (u) {
	    copy(u, u + _length, uv);
	}
	else if (upper) {
	    copy(upper, upper + _length, uv);
	}
    }
    _dist->randomSample(_data + _length * chain, _length,
			_parameters[chain], _dims, lv, uv, rng);

    delete [] lv;
    delete [] uv;
}  

bool ArrayStochasticNode::checkParentValues(unsigned int chain) const
{
    return _dist->checkParameterValue(_parameters[chain], _dims);
}

    /*
StochasticNode * 
ArrayStochasticNode::clone(vector<Node const *> const &parameters,
			   Node const *lower, Node const *upper) const
{
    return new ArrayStochasticNode(_dist, parameters, lower, upper);
}
    */
    
unsigned long ArrayStochasticNode::df() const
{
    return _dist->df(_dims);
}

void ArrayStochasticNode::sp(double *lower, double *upper, unsigned long length,
			     unsigned int chain) const
{
    _dist->support(lower, upper, length, _parameters[chain], _dims);
}

    double ArrayStochasticNode::KL(unsigned int ch1, unsigned int ch2,
				   RNG *rng, unsigned int nrep) const
    {
	if (lowerBound() && !lowerBound()->isFixed()) {
	    return JAGS_POSINF;
	}
	if (upperBound() && !upperBound()->isFixed()) {
	    return JAGS_POSINF;
	}
	return _dist->KL(_parameters[ch1], _parameters[ch2], _dims,
			 lowerLimit(ch1), upperLimit(ch1),
			 rng, nrep);
    }
    
} //namespace jags
