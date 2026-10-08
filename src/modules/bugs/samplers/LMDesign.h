#ifndef LM_DESIGN_H_
#define LM_DESIGN_H_

#include "SingletonForwardDiff.h"

namespace jags {
    namespace bugs {
	
	/**
	 * @short Calculate LM design matrix by automatic differentiation
	 *
	 * This is a wrapper around SingletonForward diff with a
	 * simplified constructor. It is used to calculate the columns
	 * of the design matrix for a linear model.
	 */
	class LMDesign : public SingletonForwardDiff
	{
	public:
	    /**
	     * Constructor.
	     *
	     * The dependent nodes passed to the SingletonForwardDiff
	     * constructor are the linear predictors - the first
	     * parent of each stochastic child of the argument gv.
	     *
	     * @param gv Singleton graph view representing a linear model.
	     * The constructor checks that the deterministic children of
	     * gv are all linear functions of the sampled node and that the
	     * stochastic children only depend on the sampled node through
	     * the first parent.
	     *
	     */
	    LMDesign(SingletonGraphView const *gv);
	};
	
    } /* namespace bugs */
} /* namespace jags */

#endif /* LM_DESIGN_H_ */
