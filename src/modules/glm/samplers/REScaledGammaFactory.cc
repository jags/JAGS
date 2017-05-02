#include <config.h>

#include "REScaledGammaFactory.h"
#include "REScaledGamma.h"

#include <graph/StochasticNode.h>

using std::vector;

namespace jags {
    namespace glm {
    
	REScaledGammaFactory::REScaledGammaFactory()
	    : REFactory("glm::REScaledGamma")
	{}

	REScaledGammaFactory::~REScaledGammaFactory()
	{}

	bool REScaledGammaFactory::canSample(StochasticNode *snode) const
	{
	    return snode->distribution()->name() == "dscaled.gamma";
	}

	REMethod *
	REScaledGammaFactory::newMethod(
	    SingletonGraphView const *tau,
	    GraphView const *eps,
	    vector<SingletonGraphView const *> const & sub_eps,
	    vector<Outcome *> const &outcomes,
	    unsigned int chain) const
	{
	    return new REScaledGamma(tau, eps, sub_eps, outcomes, chain);
	}

    } // namespace glm
} //namespace jags

