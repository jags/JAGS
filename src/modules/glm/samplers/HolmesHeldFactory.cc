#include <config.h>

#include "BinaryProbit.h"
#include "BinaryLogit.h"
#include "OrderedLogit.h"
#include "OrderedProbit.h"

#include "HolmesHeldFactory.h"
#include "HolmesHeld.h"
#include "HolmesHeldGibbs.h"

#include <module/ModuleError.h>
#include <sampler/SingletonGraphView.h>

using std::vector;

namespace jags {
namespace glm {

    HolmesHeldFactory::HolmesHeldFactory()
	: GLMFactory("glm::Holmes-Held")
    {}

    bool HolmesHeldFactory::checkOutcome(StochasticNode const *snode) const
    {
	return (BinaryProbit::canRepresent(snode) ||
		BinaryLogit::canRepresent(snode) ||
		OrderedLogit::canRepresent(snode) ||
		OrderedProbit::canRepresent(snode));
    }

    
    GLMMethod *
    HolmesHeldFactory::newMethod(GraphView const *view,
			     vector<SingletonGraphView const *> const &subviews,
			     unsigned int chain, bool gibbs) const
    {
	vector<Outcome*> outcomes;

	vector<StochasticNode *>::const_iterator p;
	for (p = view->stochasticChildren().begin();
	     p != view->stochasticChildren().end(); ++p)
	{
	    Outcome *outcome = 0;
	    if (BinaryProbit::canRepresent(*p)) {
		outcome = new BinaryProbit(*p, chain);
	    }
	    else if (BinaryLogit::canRepresent(*p)) {
		outcome = new BinaryLogit(*p, chain);
	    }
	    else if (OrderedLogit::canRepresent(*p)) {
		outcome = new OrderedLogit(*p, chain);
	    }
	    else if (OrderedProbit::canRepresent(*p)) {
		outcome = new OrderedProbit(*p, chain);
	    }

	    else {
		throwLogicError("Invalid outcome in HolmesHeldFactory");
	    }
	    outcomes.push_back(outcome);
	}

	if (gibbs) {
	    return new HolmesHeldGibbs(view, subviews, outcomes, chain);
	}
	else {
	    return new HolmesHeld(view, subviews, outcomes, chain);
	}
	
    }

}}
