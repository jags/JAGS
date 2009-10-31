#include <config.h>

#include <string>

#include "ProbitFactory.h"
#include "Probit.h"

#include <graph/StochasticNode.h>
#include <graph/LinkNode.h>
#include <distribution/Distribution.h>
#include <function/InverseLinkFunc.h>

using std::string;
using std::vector;

namespace glm {

    ProbitFactory::ProbitFactory()
	: GLMFactory("Probit")
    {}

    bool ProbitFactory::checkOutcome(StochasticNode const *snode,
				     LinkNode const *lnode) const
    {
	if (!lnode || lnode->link()->linkName() != "probit")
	    return false;

	Node const *N = 0;

	switch(GLMMethod::getFamily(snode)) {
	case GLM_BERNOULLI:
	    return true;
	case GLM_BINOMIAL:
	    N = snode->parents()[1];
	    if (N->length() != 1)
		return false;
	    if (!N->isObserved())
		return false;
	    if (N->value(0)[0] != 1)
		return false;
	    return true;
	default:
	    return false;
	}
    }
    
    GLMMethod *
    ProbitFactory::newMethod(Updater const *updater,
			     vector<Updater const *> const &sub_updaters,
			     unsigned int chain) const
    {
	return new Probit(updater, sub_updaters, chain);
    }

}
