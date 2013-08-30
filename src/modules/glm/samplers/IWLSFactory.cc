#include <config.h>

#include <graph/StochasticNode.h>

#include "IWLSFactory.h"
#include "IWLS.h"
#include <graph/LinkNode.h>

using std::vector;
using std::string;

namespace glm {

    IWLSFactory::IWLSFactory()
	: GLMFactory("glm::IWLS")
    {}

    bool IWLSFactory::checkOutcome(StochasticNode const *snode,
				   LinkNode const *lnode) const
    {
	GLMFamily family = GLMMethod::getFamily(snode);
	if (family == GLM_NORMAL) {
	    return lnode == 0;
	}
	else if (!lnode) {
	    return false;
	}
	else {
	    string link = lnode->linkName();
	    switch(family) {
	    case GLM_BERNOULLI: case GLM_BINOMIAL:
		return link == "probit" || link == "logit";
	    case GLM_POISSON:
		return link == "log";
	    case GLM_UNKNOWN: case GLM_NORMAL:
		return false;
	    }
	}
	return false; //icpc -Wall
    }
    
    GLMMethod *
    IWLSFactory::newMethod(GraphView const *view,
			     vector<GraphView const *> const &sub_views,
			     unsigned int chain) const
    {
	return new IWLS(view, sub_views, chain);
    }
    
    bool IWLSFactory::canSample(StochasticNode const *snode) const
    {
	vector<Node const *> const &parents = snode->parents();
	for (unsigned int i = 0; i < parents.size(); ++i) {
	    if (!parents[i]->isObserved())
		return false;
	}
	return !isBounded(snode);
    }

    bool IWLSFactory::fixedDesign() const
    {
	return true;
    }
}
