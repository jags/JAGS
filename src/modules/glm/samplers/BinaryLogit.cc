#include "BinaryLogit.h"
#include "Classify.h"

#include "KS.h"

#include <graph/StochasticNode.h>
#include <rng/TruncatedNormal.h>
#include <rng/RNG.h>

#include <cmath>

using std::exp;
using std::log;
using std::sqrt;

#define REG_PENALTY 0.001

//Left truncated logistic distribution
static double llogit(double left, jags::RNG *rng, double mu)
{
    double qleft = 1/(1 + exp(mu-left));
    double x = qleft + (1 - qleft) * rng->uniform();
    return mu + log(x) - log(1 - x);
}

//Right truncated logistic distribution
static double rlogit(double right, jags::RNG *rng, double mu)
{
    double qright = 1/(1 + exp(mu-right));
    double x = qright * rng->uniform();
    return mu + log(x) - log(1 - x);
}


namespace jags {
namespace glm {

    BinaryLogit::BinaryLogit(StochasticNode const *snode, unsigned int chain)
	: Outcome(snode, chain), _y(snode->value(chain)[0]), _z(0), _tau(1), _sigma2(1)
    {
	//fixme: sanity checks on snode
    }

    double BinaryLogit::value() const 
    {
	return _z;
    }

    double BinaryLogit::precision() const 
    {
	return _tau;
    }
    
    void BinaryLogit::update(RNG *rng)
    {
	/* Albert Chib update */
	if (_y) {
	    _z = llogit(0, rng, _lp);
	}
	else {
	    _z = rlogit(0, rng, _lp);
	}
	_sigma2 = sample_lambda(_z - _lp, rng);
	_tau = REG_PENALTY + 1/_sigma2;

    }

    void BinaryLogit::update(double mean, double var, RNG *rng)
    {
	/* Holmes-Held update */
	if (_y) {
	    _z = lnormal(0, rng, mean, sqrt(var + _sigma2));
	}
	else {
	    _z = rnormal(0, rng, mean, sqrt(var + _sigma2));
	}
    }

    bool BinaryLogit::canRepresent(StochasticNode const *snode)
    {
	//FIXME: Big overlap with BinaryProbit::canRepresent

	Node const *N = 0;

	switch(getFamily(snode)) {
	case GLM_BERNOULLI:
            break;
	case GLM_BINOMIAL:
            // We can also model binomial distribution if the denominator is fixed to be 1.                                                                                    
            N = snode->parents()[1];
            if (N->length() != 1)
		return false;
            if (!N->isFixed())
		return false;
            if (N->value(0)[0] != 1)
		return false;
            break;
	default:
            return false;
	}

	return getLink(snode) == LNK_LOGIT;
    }

}}
