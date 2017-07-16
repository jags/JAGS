#include <config.h>

#include "REGammaFactory2.h"
#include "REGamma2.h"

#include <graph/StochasticNode.h>

using std::vector;

namespace jags {
    namespace glm {
    
	REGammaFactory2::REGammaFactory2()
	    : REFactory2("glm::REGamma2")
	{}

	REGammaFactory2::~REGammaFactory2()
	{}

	bool REGammaFactory2::canSample(StochasticNode *snode) const
	{
	    return snode->distribution()->name() == "dgamma";
	}

	REMethod2 *
	REGammaFactory2::newMethod(SingletonGraphView const *tau,
				   GLMMethod const *glmmethod) const
	{
	    return new REGamma2(tau, glmmethod);
	}

    } // namespace glm
} //namespace jags

