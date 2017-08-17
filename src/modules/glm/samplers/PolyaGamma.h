#ifndef POLYA_GAMMA_H_
#define POLYA_GAMMA_H_

#include <config.h>

#include "Outcome.h"
#include <graph/StochasticNode.h>

namespace jags {

    class StochasticNode;
    
    namespace glm {
    
	/*
	 * @short Binary outcome with logistic link
	 *
	 * Outcome for binary logistic regression models, based on 
	 * FIXME
	 *
	 */
	class PolyaGamma : public Outcome
	{
	    double const &_y;
	    double const &_n;
	    double _tau;
	  public:
	    PolyaGamma(StochasticNode const *snode, unsigned int chain);
	    double value() const;
	    double precision() const;
	    void update(RNG *rng);
	    static bool canRepresent(StochasticNode const *snode);
	};

    }
}

#endif /* POLYA_GAMMA_H_ */
