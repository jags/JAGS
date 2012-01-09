#include <config.h>

#include "HolmesHeldBlockFactory.h"
#include "HolmesHeldBlock.h"

using std::vector;

namespace glm {

    HolmesHeldBlockFactory::HolmesHeldBlockFactory()
	: BinaryFactory("glm::Holmes-Held-Block", false)
    {}

    BinaryGLM *
    HolmesHeldBlockFactory::newBinary(GraphView const *view,
				 vector<GraphView const *> const &sub_views,
				 unsigned int chain) const
    {
	return new HolmesHeldBlock(view, sub_views, chain);
    }

    
    bool HolmesHeldBlockFactory::fixedOutcome() const
    {
	return false;
    }

    bool HolmesHeldBlockFactory::fixedDesign() const
    {
	return false;
    }
}
