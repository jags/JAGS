#include <config.h>

#include "LogisticLinear.h"
#include "Classify.h"
#include "KS.h"

#include <graph/StochasticNode.h>
#include <rng/RNG.h>

#include <cmath>

using std::sqrt;

#define REG_PENALTY 0.001


namespace jags {
    namespace glm {

	LogisticLinear::LogisticLinear(StochasticNode const *snode,
				       unsigned int chain)
	    : Outcome(snode, chain),
	      _value(snode->value(chain)[0]),
	      _mean(snode->parents()[1]->value(chain)[0]),
	      _precision(snode->parents()[1]->value(chain)[0]),
	      _lambda(1)
	{
	}
	
	void LogisticLinear::update(RNG * rng)
	{
	    double Z = (_value - _mean) * sqrt(_precision);
	    _lambda = sample_lambda(Z, rng);
	}

	double LogisticLinear::value() const 
	{
	    return _value;
	}

	double LogisticLinear::precision() const 
	{
	    return _precision * (REG_PENALTY + 1/_lambda);
	}

	bool LogisticLinear::canRepresent(StochasticNode const *snode)
	{
	    return getFamily(snode) == GLM_LOGISTIC &&
		getLink(snode) == LNK_LINEAR;
	}
    }
}
