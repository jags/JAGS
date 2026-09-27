#ifndef SINGLETON_FORWARD_DIFF_H_
#define SINGLETON_FORWARD_DIFF_H_

#include <sampler/SingletonGraphView.h>

namespace jags {
    namespace bugs {

	/**
	 * @short Automatic differentiation by forward accumulation
	 */
	class SingletonForwardDiff {
	private:
	    SingletonGraphView const *_gv;
	    std::vector<unsigned long> _edges1;
	    std::vector<std::vector<unsigned long>> _edges2;
	    std::vector<unsigned long> _terminal_nodes;
	    unsigned long _nrow;
	public:
	    SingletonForwardDiff(SingletonGraphView const *gv);
	    /**
	     * Calculate derivatives of the terminal deterministic
	     * nodes with respect to the sampled node. The terminal
	     * nodes are the direct parents of the stochastic children
	     * of the GraphView.
	     *
	     * @param coef Pointer to the start of an array which is
	     * assumed to be of size nrow times ncol.
	     *
	     * @param chain Chain number to use for the derivative
	     */
	    void calCoef(double *coef, unsigned int chain) const;
	    /**
	     * Number of rows of the coefficient matrix, equal to the
	     * sum of the lengths of the terminal nodes */
	    unsigned long nrow() const;
	    /**
	     * Number of columns of the coefficient matrix, equal to
	     * the length of the sampled node.
	     */
	    unsigned long ncol() const;
	};
    }
}

#endif /* SINGLETON_FORWARD_DIFF_H_ */
