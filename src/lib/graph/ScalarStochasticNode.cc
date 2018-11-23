#include <config.h>
#include <graph/ScalarStochasticNode.h>
#include <graph/NodeError.h>
#include <distribution/ScalarDist.h>
#include <util/nainf.h>

#include <vector>
#include <string>

using std::vector;
using std::string;

namespace jags {

ScalarStochasticNode::ScalarStochasticNode(ScalarDist const *dist,
					   unsigned int nchain,
					   vector<Node const *> const &params,
					   Node const *lower, Node const *upper)
    : StochasticNode(vector<unsigned long>(1,1), nchain, dist, params, lower, upper),
      _dist(dist)
{
    for(vector<Node const *>::const_iterator p = params.begin();
	p != params.end(); ++p)
    {
	if ((*p)->length()==0) {
	    string msg("Invalid zero-length parameter in distribution ");
	    msg.append(dist->name());
	    throw NodeError(*p, msg);
	}
	else if ((*p)->length() > 1) {
	    string msg("Invalid vector parameter in distribution ");
	    msg.append(dist->name());
	    throw NodeError(*p, msg);
	}
    }
}

double ScalarStochasticNode::logDensity(unsigned int chain, PDFType type) const
{
    if(!_dist->checkParameterValue(_parameters[chain]))
	return JAGS_NEGINF;

    double const *l = lowerLimit(chain);
    double const *u = upperLimit(chain);
    if (l && u && *l > *u) return JAGS_NEGINF;
    
    return _dist->logDensity(_data[chain], type, _parameters[chain], l, u);
}

void ScalarStochasticNode::randomSample(RNG *rng, unsigned int chain)
{
    double const *l = lowerLimit(chain);
    double const *u = upperLimit(chain);
    if (l && u && *l > *u) throw NodeError(this, "Inconsistent bounds");

    _data[chain] = _dist->randomSample(_parameters[chain], l, u, rng);
}  

void ScalarStochasticNode::truncatedSample(RNG *rng, unsigned int chain,
					   double const *lb, double const *ub)
{
    double const *l = lowerLimit(chain);
    if (lb) {
	if (l == nullptr || *lb < *l)
	    l = lb;
    }
    double const *u = upperLimit(chain);
    if (ub) {
	if (u == nullptr || *ub > *u)
	    u = ub;
    }
    if (l && u && *l > *u) throw NodeError(this, "Inconsistent bounds");
    
    _data[chain] = _dist->randomSample(_parameters[chain], l, u, rng);
}  

bool ScalarStochasticNode::checkParentValues(unsigned int chain) const
{
    double const *l = lowerLimit(chain);
    double const *u = upperLimit(chain);
    if (l && u && *l > *u) {
	return false; //Inconsistent bounds
    }
    else {
	return _dist->checkParameterValue(_parameters[chain]);
    }
}

void ScalarStochasticNode::sp(double *lower, double *upper, unsigned long,
			      unsigned int chain) const
{
    *lower = _dist->l(_parameters[chain]);
    *upper = _dist->u(_parameters[chain]);
}
    
    /*
StochasticNode * 
ScalarStochasticNode::clone(vector<Node const *> const &parameters,
			    Node const *lower, Node const *upper) const
{
    return new ScalarStochasticNode(_dist, parameters, lower, upper);
}
    */
    
    double ScalarStochasticNode::KL(unsigned int ch1, unsigned int ch2,
				    RNG *rng, unsigned int nrep) const
    {
	if (lowerBound() || upperBound()) {
	    Node const *ll = lowerBound();
	    Node const *uu = upperBound();
	    if (ll && !ll->isFixed()) {
		return JAGS_POSINF;
	    }
	    if (uu &&  !uu->isFixed()) {
		return JAGS_POSINF;
	    }
	    return _dist->KL(_parameters[ch1], _parameters[ch2],
			     lowerLimit(ch1), upperLimit(ch1),
			     rng, nrep);
	}
	else {
	    double kl =  _dist->KL(_parameters[ch1], _parameters[ch2]);
	    if (kl == JAGS_NA) {
		return _dist->KL(_parameters[ch1], _parameters[ch2],
				 nullptr, nullptr, rng, nrep);
	    }
	    else {
		return kl;
	    }
	}
    }

} //namespace jags
