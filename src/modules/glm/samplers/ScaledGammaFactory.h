#ifndef SCALED_GAMMA_FACTORY_H_
#define SCALED_GAMMA_FACTORY_H_

#include <sampler/SingletonFactory.h>

namespace jags {
    namespace glm {

	/**
	 * @short Factory for scaled precision parameters 
	 *
	 * @see ScaledGamma
	 */
	class ScaledGammaFactory : public SingletonFactory
	{
	  public:
	    ~ScaledGammaFactory();
	    bool canSample(StochasticNode *snode, Graph const &graph) const;
	    Sampler *makeSampler(StochasticNode *snode, Graph const &g) const;
	    std::string name() const;
	};
	
    }
}

#endif /* SCALED_GAMMA_FACTORY_H_ */
