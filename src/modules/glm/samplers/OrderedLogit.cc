#include "OrderedLogit.h"
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

//Interval truncated logistic distribution
static double ilogit(double left, double right, jags::RNG *rng, double mu)
{
    double qleft = 1/(1 + exp(mu-left));
    double qright = 1/(1 + exp(mu-right));
    double x = qleft + (qright - qleft) * rng->uniform();
    return mu + log(x) - log(1 - x);
}


namespace jags {
namespace glm {

    OrderedLogit::OrderedLogit(StochasticNode const *snode, unsigned int chain)
	: Outcome(snode, chain), _y(snode->value(chain)[0]),
	  _cuts(snode->parents()[1]->value(chain)),
	  _ncut(snode->parents()[1]->length()),
	  _z(0), _tau(1), _sigma2(1)
    {
	//fixme: sanity checks on snode
    }

    double OrderedLogit::value() const 
    {
	return _z;
    }

    double OrderedLogit::precision() const 
    {
	return _tau;
    }
    
    void OrderedLogit::update(RNG *rng)
    {
	/* Albert Chib update */
	int x = static_cast<int>(_y) - 1;
	if (x == 0) {
	    _z = rlogit(_cuts[0], rng, _lp);
	}
	else if (x == _ncut) {
	    _z = llogit(_cuts[_ncut - 1], rng, _lp);
	}
	else {
	    _z = ilogit(_cuts[x-1], _cuts[x], rng, _lp);
	}
	
	_sigma2 = sample_lambda(_z - _lp, rng);
	_tau = REG_PENALTY + 1/_sigma2;

    }

    void OrderedLogit::update(double mean, double var, RNG *rng)
    {
	/* Holmes-Held update */
	int x = static_cast<int>(_y) - 1;
	if (x == 0) {
	    _z = rnormal(_cuts[0], rng, mean, sqrt(var + _sigma2));
	}
	else if (x == _ncut) {
	    _z = lnormal(_cuts[_ncut - 1], rng, mean, sqrt(var + _sigma2));
	}
	else {
	    _z = inormal(_cuts[x-1], _cuts[x], rng, mean, sqrt(var + _sigma2));
	}
    }

    bool OrderedLogit::canRepresent(StochasticNode const *snode)
    {
	return getFamily(snode) == GLM_ORDLOGIT &&
	    getLink(snode) == LNK_LINEAR;
    }

}}
