#ifndef RE_SCALED_WISHART_H_
#define RE_SCALED_WISHART_H_

#include "REMethod.h"

namespace jags {
    namespace glm {

	/**
	 * @short Random effects sampler for scaled Wishart precision
	 */
	class REScaledWishart : public REMethod {
	    std::vector<double> _sigma;
	  public:
	    REScaledWishart(SingletonGraphView const *tau,
			    GraphView const *eps, 
			    std::vector<SingletonGraphView const *> const &veps,
			    std::vector<Outcome *> const &outcomes,
			    unsigned int chain);
	    void updateTau(RNG *rng);
	    void updateSigma(RNG *rng);
	};

    }
}

#endif /* RE_SCALED_WISHART_H_ */
