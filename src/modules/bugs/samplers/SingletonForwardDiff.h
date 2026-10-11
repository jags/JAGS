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
	    std::vector<unsigned long> _dependent_nodes, _gindex;
	    unsigned long _nrow;
	public:
	    /**
	     * Constructor.
	     *
	     * @param gv Singleton graph view. The derivative will be
	     * calculated with respect to the sampled node.
	     *
	     * @param dependent_nodes Vector of nodes representing
	     * dependent variables. The derivative of each terminal
	     * node is calculated with respect to the sampled node.
	     */
	    SingletonForwardDiff(SingletonGraphView const *gv,
				 std::vector<Node const *> const &dependent_nodes);
	    /**
	     * Calculate derivatives of the dependent nodes with
	     * respect to the sampled node.
	     *
	     * @param coef Pointer to the start of an array which is
	     * assumed to be of size nrow times ncol.
	     *
	     * @param chain Chain number to use for the derivative
	     *
	     * @param transpose Boolean flag. If true then the transpose of the
	     * gradient matrix is returned.
	     */
	    void calCoef(double *coef, unsigned int chain, bool transpose=false) const;
	    /**
	     * Number of rows of the coefficient matrix, equal to the
	     * sum of the lengths of dependent nodes
	     */
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
