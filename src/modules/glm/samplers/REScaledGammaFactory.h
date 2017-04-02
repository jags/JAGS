#ifndef RE_SCALED_GAMMA_FACTORY_H_
#define RE_SCALED_GAMMA_FACTORY_H_

#include "REFactory.h"

namespace jags {
    namespace glm {

	/**
	 * @short Factory for random effects with scaled gamma precision
	 *
	 * @see REScaledGamma
	 */
	class REScaledGammaFactory : public REFactory
	{
	  public:
	    REScaledGammaFactory();
	    ~REScaledGammaFactory();
	    bool canSample(StochasticNode *snode) const;
	    REMethod * newMethod(SingletonGraphView const *tau,
				 GraphView const *eps,
				 std::vector<SingletonGraphView const *> const & veps,
				 std::vector<Outcome*> const &outcomes,
				 unsigned int chain) const;
	};
	
    }
}

#endif /* RE_SCALED_GAMMA_FACTORY_H_ */
