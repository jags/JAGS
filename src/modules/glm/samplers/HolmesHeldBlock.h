#ifndef HOLMES_HELD_BLOCK_H_
#define HOLMES_HELD_BLOCK_H_

#include "BinaryGLM.h"

namespace glm {

    /**
     * @short Holmes Held sampler for binary GLMs
     *
     * Sampler for probit and logistic regression models with binary
     * outcome data.
     */
    class HolmesHeldBlock : public BinaryGLM {
      public:
	/**
	 * Constructor.
	 * 
	 * @see GLMMethod#GLMMethod
	 */
	HolmesHeldBlock(GraphView const *view, 
		      std::vector<GraphView const *> const &sub_views,
		      unsigned int chain);
	/**
	 * Does block updating of latent auxiliary variables and their
	 * variances
	 */
	void updateAuxiliary(cholmod_dense *u1, cholmod_factor *F, RNG *rng);
	/**
	 * Returns "Holmes-Held Block"
	 */
	std::string name() const;
	/**
	 * Joint update of model parameters and auxiliary variables
	 */
	void update(RNG *rng);
    };
    
}

#endif /* HOLMES_HELD_BLOCK_H_ */
