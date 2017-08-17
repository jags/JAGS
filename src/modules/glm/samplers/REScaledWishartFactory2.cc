#include <config.h>

#include "REScaledWishartFactory2.h"
#include "REScaledWishart2.h"

#include <graph/StochasticNode.h>

using std::vector;

namespace jags {
    namespace glm {
    
	REScaledWishartFactory2::REScaledWishartFactory2()
	    : REFactory2("glm::REScaledWishart2")
	{}

	REScaledWishartFactory2::~REScaledWishartFactory2()
	{}

	bool REScaledWishartFactory2::canSample(StochasticNode *snode) const
	{
	    return snode->distribution()->name() == "dscaled.wishart";
	}

	REMethod2 *
	REScaledWishartFactory2::newMethod(SingletonGraphView const *tau,
					   GLMMethod const *glmmethod) const
	{
	    return new REScaledWishart2(tau, glmmethod);
	}

    } // namespace glm
} //namespace jags

