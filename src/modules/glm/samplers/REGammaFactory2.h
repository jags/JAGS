#ifndef RE_GAMMA_FACTORY2_H_
#define RE_GAMMA_FACTORY2_H_

#include "REFactory2.h"

namespace jags {
    namespace glm {

	/**
	 * @short Factory for random effects with gamma precision
	 *
	 * @see REGamma2
	 */
	class REGammaFactory2 : public REFactory2
	{
	  public:
	    REGammaFactory2();
	    ~REGammaFactory2();
	    bool canSample(StochasticNode *snode) const;
	    REMethod2 *	newMethod(SingletonGraphView const *tau,
				  GLMMethod const *randef) const;
	};
	
    }
}

#endif /* RE_GAMMA_FACTORY_H_ */
