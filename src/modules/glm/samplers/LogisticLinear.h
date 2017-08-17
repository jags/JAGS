#ifndef LOGISTIC_LINEAR_H_
#define LOGISTIC_LINEAR_H_

#include <config.h>
#include "Outcome.h"

#include <graph/StochasticNode.h>

namespace jags {

    class StochasticNode;
    
    namespace glm {

	class LogisticLinear : public Outcome
	{
	    double const &_value;
	    double const &_mean;
	    double const &_precision;
	    double _lambda;
	  public:
	    LogisticLinear(StochasticNode const *snode, unsigned int chain);
	    double value() const;
	    double precision() const;
	    void update(RNG *rng);
	    static bool canRepresent(StochasticNode const *snode);
	};

    }
}

#endif /* LOGISTIC_LINEAR_H_ */
