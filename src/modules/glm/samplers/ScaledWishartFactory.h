#ifndef SCALED_WISHART_FACTORY_H_
#define SCALED_WISHART_FACTORY_H_

#include <sampler/SingletonFactory.h>

namespace jags {
    namespace glm {

	/**
	 * @short Factory for scaled precision parameters 
	 *
	 * @see ScaledWishart
	 */
	class ScaledWishartFactory : public SingletonFactory
	{
	  public:
	    ~ScaledWishartFactory();
	    bool canSample(StochasticNode *snode, Graph const &graph) const;
	    Sampler *makeSampler(StochasticNode *snode, Graph const &g) const;
	    std::string name() const;
	};
	
    }
}

#endif /* SCALED_WISHART_FACTORY_H_ */
