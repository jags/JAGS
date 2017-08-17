#include <config.h>

#include "REScaledGammaFactory2.h"
#include "REScaledGamma2.h"

#include <graph/StochasticNode.h>

using std::vector;

namespace jags {
    namespace glm {
    
	REScaledGammaFactory2::REScaledGammaFactory2()
	    : REFactory2("glm::REScaledGamma2")
	{}

	REScaledGammaFactory2::~REScaledGammaFactory2()
	{}

	bool REScaledGammaFactory2::canSample(StochasticNode *snode) const
	{
	    return snode->distribution()->name() == "dscaled.gamma";
	}

	REMethod2 *
	REScaledGammaFactory2::newMethod(SingletonGraphView const *tau,
					 GLMMethod const *glmmethod) const
	{
	    return new REScaledGamma2(tau, glmmethod);
	}

    } // namespace glm
} //namespace jags

