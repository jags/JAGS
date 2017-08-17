#include <config.h>

#include "IWLSFactory.h"
#include "IWLS.h"
#include "NormalLinear.h"
#include "IWLSOutcome.h"

#include <graph/StochasticNode.h>
#include <sampler/GraphView.h>

using std::vector;
using std::string;

namespace jags {
namespace glm {

    IWLSFactory::IWLSFactory()
	: GLMFactory("glm::IWLS")
    {}

    bool IWLSFactory::checkOutcome(StochasticNode const *snode) const
    {
	return NormalLinear::canRepresent(snode) || 
	    IWLSOutcome::canRepresent(snode);
    }
    
    GLMMethod *
    IWLSFactory::newMethod(GraphView const *view,
			   vector<SingletonGraphView const *> const &sub_views,
			   unsigned int chain, bool gibbs) const
    {
        bool linear = true;
        vector<Outcome*> outcomes;

        for (vector<StochasticNode *>::const_iterator p = view->stochasticChildren().begin();
             p != view->stochasticChildren().end(); ++p)
        {
            Outcome *outcome = 0;
            if (NormalLinear::canRepresent(*p)) {
                outcome = new NormalLinear(*p, chain);
            }
            else if (IWLSOutcome::canRepresent(*p)) {
                outcome = new IWLSOutcome(*p, chain);
                linear = false;
            }
            outcomes.push_back(outcome);
        }

        if (linear) {
            return new GLMBlock(view, sub_views, outcomes, chain);
        }
	
	return new IWLS(view, sub_views, outcomes, chain);
    }
    
    bool IWLSFactory::canSample(StochasticNode const *snode) const
    {
	vector<Node const *> const &parents = snode->parents();
	for (unsigned int i = 0; i < parents.size(); ++i) {
	    if (!parents[i]->isFixed())
		return false;
	}
	return !isBounded(snode);
    }

    bool IWLSFactory::fixedDesign() const
    {
	return true;
    }

}}
