#ifndef RE_SCALED_GAMMA_FACTORY2_H_
#define RE_SCALED_GAMMA_FACTORY2_H_

#include "REFactory2.h"

namespace jags {
    namespace glm {

	/**
	 * @short Factory for random effects with scaled gamma precision
	 *
	 * @see REScaledGamma2
	 */
	class REScaledGammaFactory2 : public REFactory2
	{
	  public:
	    REScaledGammaFactory2();
	    ~REScaledGammaFactory2();
	    bool canSample(StochasticNode *snode) const;
	    REMethod2 *	newMethod(SingletonGraphView const *tau,
				  GLMMethod const *randef) const;
	};
	
    }
}

#endif /* RE_SCALED_GAMMA_FACTORY2_H_ */
