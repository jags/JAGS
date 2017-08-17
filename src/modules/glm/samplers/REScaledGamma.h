#ifndef RE_SCALED_GAMMA_H_
#define RE_SCALED_GAMMA_H_

#include "REMethod.h"

namespace jags {
    namespace glm {

	/**
	 * @short Random effects sampler for scaled gamma precision
	 */
	class REScaledGamma : public REMethod {
	    double _sigma;
	  public:
	    REScaledGamma(SingletonGraphView const *tau,
			  GraphView const *eps, 
			  std::vector<SingletonGraphView const *> const &veps,
			  std::vector<Outcome *> const &outcomes,
			  unsigned int chain);
	    void updateTau(RNG *rng);
	    void updateSigma(RNG *rng);
	};

    }
}

#endif /* RE_SCALED_GAMMA_H_ */
