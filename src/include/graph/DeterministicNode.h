#ifndef DETERMINISTIC_NODE_H_
#define DETERMINISTIC_NODE_H_

#include <graph/Node.h>

#include <set>

namespace jags {

/**
 * @short Closed classes of deterministic nodes
 *
 * A deterministic node may be considered as a function of any of its
 * deterministic ancestors.  Some classes of function are closed under
 * composition: if g() and f() both belong to the class, then so does
 * g(f()). This property is used by some Samplers to determine whether
 * a node, or block of nodes, can be sampled.
 *
 * DNODE_ADDITIVE denotes a function f of the single ancestor nodes X
 * such that sum(f(X)) = a + sum(X), where sum() denotes the sum of
 * all the elements of the given node.
 *
 * DNODE_LINEAR denotes a linear function of the ancestor nodes X =
 * (X1, ... XN). A linear function takes the form A + B %*% X1 + B2
 * %*% X2 + ... * + BN %*% * XN. 
 *
 * DNODE_SCALE denotes a scale function of the single ancestor node X.
 * This is a function of the form B %*% X.
 *
 * DNODE_SCALE_MIX denotes a generalization of the scale function
 * class used in mixture models. Scale mixture functions are linear
 * functions of the form A + B %*% X, where at least one of
 * A or B is zero.
 *
 * DNODE_POWER denotes a power function of the single ancestor node X.
 * This is a function of the form Y = A * X^B.  Power
 * functions are linear functions on a log scale: i.e. if Y=A*X^B and
 * A > 0, then log(Y) is a linear function of log(X).
 *
 * A function is considered to be a fixed function if the coefficient
 * B, or coefficients B=(B1...BN) is fixed. 
 *
 * @see DeterministicNode#isClosed
 */
enum ClosedFuncClass {DNODE_ADDITIVE, DNODE_LINEAR, DNODE_SCALE,
		      DNODE_SCALE_MIX, DNODE_POWER};

/**
 * @short Base class for deterministic Nodes
 *
 * The value of a deterministic node is determined exactly by the
 * values of its parents.
 */
class DeterministicNode : public Node {
    bool _fixed;
    const std::array<int, 2> _depth;
public:
    DeterministicNode(std::vector<unsigned long> const &dim,
		      unsigned int nchain,
		      std::vector<Node const *> const &parents);
    ~DeterministicNode();
    /**
     * A DeterministicNode has the same stochastic depth as its
     * deepest parent. Its deterministic depth is one greater than its
     * deepest parent with the same stochastic depth.
     */
    std::array<int, 2> const &depth() const;
    /**
     * Random samples from a Deterministic node are not random.
     * This function simply calculates the value of the node from its
     * parent values and leaves the RNG object untouched.
     */
    void randomSample(RNG*, unsigned int nchain);
    /**
     * Deterministic nodes are not random variables. 
     */
    bool isRandomVariable() const;
    /**
     * A deterministic node is never observed.
     */
    bool isObserved(unsigned long offset) const;
    /**
     * A deterministic node is fixed if all its parents are fixed.
     */
    bool isFixed() const;
    /**
     * Tests whether the node belongs to a closed class when
     * considered as a function of a given ancestor node X, or nodes
     * X1 ... Xn.
     *
     * False negative responses are permitted: i.e. the value false
     * may be returned when the node is, in fact, in the closed class,
     * but false positives are not allowed.
     *
     * A pre-condition is that all nodes on the path from the ancestor
     * node X to the current node must be in the closed class. Thus,
     * the isClosed function should be called iteratively, and if it
     * returns false for any ancestor then you can infer that this
     * node is not in the closed class without calling the isClosed
     * function.
     *
     * @param ancestors Set containing all ancestors of the test node
     * that are in the closed class.
     *
     * @param fc Closed class to be tested.
     *
     * @param fixed When true, the closed class is restricted to the
     * sub-class of functions that are considered fixed.
     */
    virtual bool isClosed(std::set<Node const *> const &ancestors, 
			  ClosedFuncClass fc, bool fixed) const = 0;
    /*
     * Creates a copy of the deterministic node.  Supplying the parents
     * of this node as the argument creates an identical copy.
     *
     * @param parents Parents of the cloned node. 
     *
    virtual DeterministicNode * clone(std::vector<Node const *> const &parents)
	const = 0;
     */
    void unlinkParents();
    /**
     * Calculates a value for the node based on its parents' values.
     *
     * @param chain Number of chain from which to draw sample
     */
    virtual void deterministicSample(unsigned int chain) = 0;

    double logDensity(unsigned int chain, PDFType type) const;
    unsigned int df() const;
    double KL(unsigned int chain1, unsigned int chain2, RNG *rng,
	      unsigned int nrep) const;
    /**
     * Used by dumpNodeNames to gather a specific subset of node types:
     */
    inline bool isConstant() const { return false; }
    inline bool isDeterministic() const { return true; }
    inline bool isStochastic() const { return false; }
};

} /* namespace jags */

#endif /* DETERMINISTIC_NODE_H_ */
