#include <config.h>

#include "REGammaFactory.h"
#include "REGamma.h"

#include <graph/StochasticNode.h>

using std::vector;

namespace jags {
    namespace glm {
    
	REGammaFactory::REGammaFactory()
	    : REFactory("REGamma")
	{}

	REGammaFactory::~REGammaFactory()
	{}

	bool REGammaFactory::canSample(StochasticNode *snode) const
	{
	    return snode->distribution()->name() == "dgamma";
	}

	REMethod *
	REGammaFactory::newMethod(
	    SingletonGraphView const *tau,
	    GraphView const *eps,
	    vector<SingletonGraphView const *> const & sub_eps,
	    vector<Outcome *> const &outcomes,
	    unsigned int chain) const
	{
	    return new REGamma(tau, eps, sub_eps, outcomes, chain);
	}

    } // namespace glm
} //namespace jags

