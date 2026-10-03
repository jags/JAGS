#ifndef GLM_DESIGN_H_
#define GLM_DESIGN_H_

#include <sampler/SingletonGraphView.h>

namespace jags {
    namespace glm {

	/**
	 * @short Calculate GLM design matrix by automatic differentiation
	 */
	class GLMDesign {
	private:
	    SingletonGraphView const *_gv;
	    std::vector<unsigned long> _edges1;
	    std::vector<std::vector<unsigned long>> _edges2;
	    std::vector<unsigned long> _terminal_nodes;
	    unsigned long _nrow;
	public:
	    /**
	     * Constructor.
	     *
	     * @param gv Singleton graph view. The columns of the design matrix
	     * will be calculated for the sampled node.
	     *
	     * @param stoch_children Vector of stochastic children for
	     * the full design matrix. This is required to give a
	     * canonical ordering to the rows of the design matrix.
	     */
	    GLMDesign(SingletonGraphView const *gv,
		      std::vector<StochasticNode *> const &stoch_children);
	    /**
	     * Calculate the coefficients of the linear predictor which form
	     * the columns of the design matrix
	     *
	     * @param coef Pointer to the start of an array which is
	     * assumed to be of size nrow times ncol.
	     *
	     * @param chain Chain number to evaluate
	     */
	    void calCoef(double *coef, unsigned int chain) const;
	    /**
	     * Number of rows of the coefficient matrix, equal to the
	     * sum of the lengths of the stochastic children */
	    unsigned long nrow() const;
	    /**
	     * Number of columns of the coefficient matrix, equal to
	     * the length of the sampled node.
	     */
	    unsigned long ncol() const;
	};
    }
}

#endif /* GLM_DESIGN_H_ */
