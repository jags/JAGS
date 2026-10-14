#ifndef GLM_DESIGN_H_
#define GLM_DESIGN_H_

#include "SingletonForwardDiff.h"

#include <sampler/SingletonGraphView.h>

namespace jags {
    namespace glm {

	/**
	 * @short Calculate GLM design matrix by automatic differentiation
	 */
	class GLMDesign : public SingletonForwardDiff {
	public:
	    /**
	     * Constructor.
	     *
	     * @param gv Singleton graph view representing a
	     * generalized linear model.  The constructor checks that
	     * gv represents a GLM and that its stochastic children
	     * only depend on the sampled node through the first
	     * parent.
	     *
	     * @param outcomes Vector of stochastic nodes representing
	     * outcomes in the GLM. The dependent nodes passed to the
	     * SingletonForwardDiff constructor are the linear
	     * predictors for these nodes.
	     *
	     * The outcomes must be a super-set of the stochastic children
	     * of gv and must not contain any repeats.
	     */
	    GLMDesign(SingletonGraphView const *gv,
		      std::vector<StochasticNode *> const &outcomes);
	};

    }
}

#endif /* GLM_DESIGN_H_ */
