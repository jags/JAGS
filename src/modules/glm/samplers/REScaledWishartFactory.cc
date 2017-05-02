#include <config.h>

#include "REScaledWishartFactory.h"
#include "REScaledWishart.h"

#include <graph/StochasticNode.h>

using std::vector;

namespace jags {
    namespace glm {
    
	REScaledWishartFactory::REScaledWishartFactory()
	    : REFactory("glm::REScaledWishart")
	{}

	REScaledWishartFactory::~REScaledWishartFactory()
	{}

	bool REScaledWishartFactory::canSample(StochasticNode *snode) const
	{
	    return snode->distribution()->name() == "dscaled.wishart";
	}

	REMethod *
	REScaledWishartFactory::newMethod(
	    SingletonGraphView const *tau,
	    GraphView const *eps,
	    vector<SingletonGraphView const *> const & sub_eps,
	    vector<Outcome *> const &outcomes,
	    unsigned int chain) const
	{
	    return new REScaledWishart(tau, eps, sub_eps, outcomes, chain);
	}

    } // namespace glm
} //namespace jags

