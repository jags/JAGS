#ifndef HOLMES_HELD_BLOCK_FACTORY_H_
#define HOLMES_HELD_BLOCK_FACTORY_H_

#include "BinaryFactory.h"

namespace glm {

    /**
     * @short Factory object for the Holmes-Held Block sampling method
     */
    class HolmesHeldBlockFactory : public BinaryFactory
    {
    public:
	HolmesHeldBlockFactory();
	/**
	 * Returns a newly allocated object of class HolmesHeldBlock for
	 * sampling binary GLMs with probit or logistic link.
	 */
	BinaryGLM *newBinary(GraphView const *view, 
			     std::vector<GraphView const *> const &sub_views,
			     unsigned int chain) const;
	/**
	 * Returns false. 
	 */
	bool fixedOutcome() const;
	/**
	 * Returns false
	 */
	bool fixedDesign() const;
    };

}

#endif /* HOLMES_HELD_BLOCK_FACTORY_H_ */
