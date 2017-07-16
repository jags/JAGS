#ifndef RE_SCALED_WISHART2_H_
#define RE_SCALED_WISHART2_H_

#include "REMethod2.h"

namespace jags {
    namespace glm {

	/**
	 * @short Random effects sampler for scaled Wishart precision
	 */
	class REScaledWishart2 : public REMethod2 {
	    std::vector<double> _sigma;
	  public:
	    REScaledWishart2(SingletonGraphView const *tau,
			     GLMMethod const *glmmethod);
	    void updateTau(RNG *rng);
	    void updateSigma(RNG *rng);
	    bool isAdaptive() const;
	    void adaptOff();
	    bool checkAdaptation() const;
	};

    }
}

#endif /* RE_SCALED_WISHART2_H_ */
