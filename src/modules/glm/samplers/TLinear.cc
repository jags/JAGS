#include <config.h>

#include <JRmath.h>

#include "TLinear.h"
#include "Classify.h"

//#include <graph/StochasticNode.h>
#include <rng/RNG.h>

namespace jags {
    namespace glm {

	TLinear::TLinear(StochasticNode const *snode, unsigned int chain)
	    : Outcome(snode, chain),
	      _value(snode->value(chain)[0]),
	      _mean(snode->parents()[0]->value(chain)[0]),
	      _precision(snode->parents()[1]->value(chain)[0]),
	      _df(snode->parents()[2]->value(chain)[0]),
	      _lambda(1)
	{
	}
	
	void TLinear::update(RNG * rng)
	{
	    double delta = _value - _mean;
	    double r = (_df + 1.0)/2;
	    double mu = (_df + _precision * delta * delta)/2;

	    _lambda = rgamma(r, 1/mu, rng);
	}

	double TLinear::value() const 
	{
	    return _value;
	}

	double TLinear::precision() const 
	{
	    return _precision * _lambda;
	}

	bool TLinear::canRepresent(StochasticNode const *snode)
	{
	    return getFamily(snode) == GLM_T &&	getLink(snode) == LNK_LINEAR;
	}
    }
}
