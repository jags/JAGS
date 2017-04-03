#ifndef RE_GAMMA_H_
#define RE_GAMMA_H_

#include "REMethod.h"

namespace jags {
    namespace glm {

	/**
	 * @short Random effects sampler for gamma precision
	 */
	class REGamma : public REMethod {
	  public:
	    REGamma(SingletonGraphView const *tau,
		    GraphView const *eps, 
		    std::vector<SingletonGraphView const *> const &veps,
		    std::vector<Outcome *> const &outcomes,
		    unsigned int chain);
	    void updateTau(RNG *rng);
	};

    }
}

#endif /* RE_GAMMA_H_ */
