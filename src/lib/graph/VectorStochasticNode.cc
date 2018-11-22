#include <config.h>
#include <graph/VectorStochasticNode.h>
#include <graph/NodeError.h>
#include <distribution/DistError.h>
#include <distribution/VectorDist.h>
#include <util/nainf.h>
#include <util/dim.h>

#include <vector>
#include <string>
#include <algorithm>

using std::vector;
using std::string;
using std::max;
using std::min;
using std::copy;

namespace jags {

static unsigned long mkLength(VectorDist const *dist, 
			     vector<Node const *> const &parents)
{
    /* 
       Calculates length of stochastic node as a function of its
       parents
    */

    if (!checkNPar(dist, parents.size())) {
	throw DistError(dist, "Incorrect number of parameters");
    }
    vector<unsigned long> parameter_lengths(parents.size());
    for (unsigned long j = 0; j < parents.size(); ++j) {
	parameter_lengths[j] = parents[j]->length();
    }
    if (!dist->checkParameterLength(parameter_lengths)) {
	throw DistError(dist, "Invalid parameter lengths");
    }
    return dist->length(parameter_lengths);
}

static vector<unsigned long> const &
mkParameterLengths(vector<Node const *> const &parameters) {
    vector<unsigned long>  lengths(parameters.size());
    for (unsigned long j = 0; j < parameters.size(); ++j) {
        lengths[j] = parameters[j]->length();
    }
    return getUnique(lengths);
}

VectorStochasticNode::VectorStochasticNode(VectorDist const *dist,
					   unsigned int nchain,
					   vector<Node const *> const &params)
    : StochasticNode(vector<unsigned long>(1,mkLength(dist, params)), 
		     nchain, dist, params, nullptr, nullptr),
      _dist(dist), _lengths(mkParameterLengths(params))
{
    if (!dist->checkParameterLength(_lengths)) {
	throw DistError(dist, "Invalid parameter lengths");
    }
}

double VectorStochasticNode::logDensity(unsigned int chain, PDFType type) const
{
    if(!_dist->checkParameterValue(_parameters[chain], _lengths))
	return JAGS_NEGINF;
    
    return _dist->logDensity(_data + _length * chain, type,
			     _parameters[chain], _lengths);
}

void VectorStochasticNode::randomSample(RNG *rng, unsigned int chain)
{
    _dist->randomSample(_data + _length * chain, 
			_parameters[chain], _lengths, rng);
}

bool VectorStochasticNode::checkParentValues(unsigned int chain) const
{
    return _dist->checkParameterValue(_parameters[chain], _lengths);
}

    /*
StochasticNode * 
VectorStochasticNode::clone(vector<Node const *> const &parameters,
			    Node const *lower, Node const *upper) const
{
    return new VectorStochasticNode(_dist, parameters, lower, upper);
}
    */
    
unsigned long VectorStochasticNode::df() const
{
    return _dist->df(_lengths);
}

void VectorStochasticNode::sp(double *lower, double *upper,
			      unsigned long length, unsigned int chain) const
{
    _dist->support(lower, upper, _parameters[chain], _lengths);
}

    
    double VectorStochasticNode::KL(unsigned int ch1, unsigned int ch2,
				    RNG *rng, unsigned int nrep) const
    {
	double kl =  _dist->KL(_parameters[ch1], _parameters[ch2],
			       _lengths);
	if (kl == JAGS_NA) {
	    return _dist->KL(_parameters[ch1], _parameters[ch2], _lengths,
			     rng, nrep);
	}
	else {
	    return kl;
	}
    }
    
} //namespace jags
