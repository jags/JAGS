#ifndef GLM_BLOCK_H_
#define GLM_BLOCK_H_

#include "GLMMethod.h"

#include <vector>

namespace jags {

    struct RNG;
    class GraphView;
    class SingletonGraphView;

namespace glm {

    class Outcome;

    /**
     * @short Block sampler for generalized linear models.
     */
    class GLMBlock : public GLMMethod {
    public:
	/**
	 * Constructor.
	 *
	 * @param view Pointer to a GraphView object for all sampled nodes.
	 *
	 * @param sub_views Vector of pointers to SingletonGraphView
	 * objects with length equal to the number of sampled
	 * nodes. Each sub-view corresponds to a single sampled node.
	 *
	 * @param outcomes Vector of pointers to Outcome objects with length
	 * equal to the number of stochastic children of the sampled nodes.
	 * The GLMMethod objects takes ownership of each Outcome in the vector
	 * and frees the memory in the destructor.
	 * 
	 * @param chain Number of the chain (starting from 0) to which
	 * the sampling method will be applied.
	 */
	GLMBlock(GraphView const *view, 
		 std::vector<SingletonGraphView const *> const &sub_views,
		 std::vector<Outcome *> const &outcomes,
		 unsigned int chain);
	/**
	 * Updates the regression parameters by treating the GLM as a
	 * linear model (LM).  All regression parameters are updated
	 * together in a block.
	 *
	 * @param rng Random number generator used for sampling
	 */
	void update(RNG *rng);
	/**
	 * Updates auxiliary variables.  The default method does
	 * nothing.  Sampling methods that use auxiliary variables to
	 * reduce the GLM to a linear model must provide their own
	 * implementation.
	 *
	 * This function is called by updateLM.  Internally, updateLM
	 * calculates the posterior mean "mu" by solving the equation
	 * (A %*% mu = b) where "A" is the posterior precision. The
	 * same dense vector is used to hold "b" (before solving the
	 * equation) and "mu" (after solving the equation).  If
	 * updateLM is called with parameter "chol" set to false then
	 * updateAuxiliary is called before solving the equation: thus
	 * the first argument (y) should contain the canonical
	 * parameter (b).  If updateLM is called with "chol" set to
	 * true then updateAuxiliary is called after solving the
	 * equation: thus the first argument (y) should contain the
	 * posterior mean.
	 *
	 * IMPORTANT NOTE: GLMMethod uses a parameterization in which
	 * the current value of the parameters is considered the
	 * origin. The value of "y" (mu or b) may need to be adjusted
	 * for this centring by an implementation of updateAuxiliary.
	 * 
	 * @param y Dense vector which may be either the posterior
	 * mean "mu" or (A %*% mu), where "A" is the posterior
	 * precision.
	 *
	 * @param N Cholesky factorization of the posterior precision "A".
	 *
	 * @param rng Random number generator used for sampling.
	 */
	virtual void updateAuxiliary(cholmod_dense *y, cholmod_factor *N, 
				     RNG *rng);
    };

}}

#endif /* GLM_BLOCK_H_ */
