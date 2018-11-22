#ifndef NODE_H_
#define NODE_H_

#include <list>
#include <vector>
#include <string>
#include <array>

#include <distribution/Distribution.h>
// Required for PDFtype enum

namespace jags {

struct RNG;
class StochasticNode;
class DeterministicNode;
class Graph;

/**
 * @short Node in a directed acyclic graph
 */
class Node {
    std::vector<Node const *> _parents;
    std::list<StochasticNode*> *_stoch_children;
    std::list<DeterministicNode *> *_dtrm_children;

    /* Forbid copying of Node objects */
    Node(Node const &orig);
    Node &operator=(Node const &rhs);

protected:
    std::vector<unsigned long> const &_dim;
    const unsigned long _length;
    const unsigned int _nchain;
    double *_data;

public:
    /**
     * Constucts a Node with no parents.
     * @param dim Dimension of new Node.
     * @param nchain Number of chains that the Node will contain data on.
     */
    Node(std::vector<unsigned long> const &dim, unsigned int nchain);
    /**
     * Constructs a node with parents.  Each parent must contain data on
     * the same number of chains. Subclasses of Node may give specific
     * meaning to the ordering of the parents.
     *
     * @param dim Dimension of new Node.
     *
     * @param parents vector of parent nodes. A node may not be its own
     * parent.
     */
    Node(std::vector<unsigned long> const &dim, unsigned int nchain,
	 std::vector<Node const *> const &parents);
    /**
     * Destructor. 
     */
    virtual ~Node();
    /**
     * Number of chains.
     */ 
    unsigned int nchain() const;
    /**
     * Vector of parents.
     */
    std::vector<Node const *> const &parents() const;
    /**
     * Returns the topological depth of a node as an array of length
     * 2. The first element is the stochastic depth, the second is the
     * deterministic depth.
     */
    virtual std::array<int, 2> const &depth() const = 0;
    /**
     * Draws a random sample from the node's prior distribution.
     * @param rng Pointer to random number generator
     * @param chain Number of chain from which to draw sample
     */
    virtual void randomSample(RNG *rng, unsigned int chain) = 0;
    /**
     * Checks whether the parents of the Node have valid values.
     */
    virtual bool checkParentValues(unsigned int chain) const = 0;
    /**
     * Returns the stochastic children of the node
     */
    std::list<StochasticNode*> const *stochasticChildren();
    /**
     * Returns the deterministic children of the node
     */
    std::list<DeterministicNode*> const *deterministicChildren();
    /**
     * Initializes the node for the given chain. The value array of a
     * newly constructed Node consists of missing values (denoted by
     * the special value JAGS_NA).  This function sets the value of
     * the node by forward sampling from its parents.  If the Node has
     * previously had its value set, the function will do nothing and
     * return the value true.  Initialization will fail if any of the
     * parent nodes is uninitialized, and in this case the return
     * value is false.
     *
     * @param rng RNG to use for random sampling
     * 
     * @param chain Index number of chain to initialize.
     *
     * @returns a logical value indicating success
     */
    bool initialize(RNG *rng, unsigned int chain);
    /**
     * Returns the BUGS-language representation of the node, based on the 
     * names of its parents
     *
     * @param parents Vector of names of parent nodes
     */
    virtual std::string 
	deparse(std::vector<std::string> const &parents) const = 0;
    /**
     * Indicates whether the Node represents a random variable.
     */
    virtual bool isRandomVariable() const = 0;
    /**
     * Indicates whether the value of the node is fixed.
     */
    virtual bool isFixed() const = 0;
    /**
     * Indicates whether an element of the node is observed.  Only
     * random variables can be observed. A node that is fixed but not
     * a random variable is not considered to be observed. A random
     * variable may be partially observed, so it is possible for this
     * function to return different values for different elements of
     * the same node.
     *
     * @param i index of element to be tested. If this exceeds the
     * length of the node then the return value is false.
     */
    virtual bool isObserved(unsigned long i) const = 0;
    /**
     * Sets the value of the node for a given chain
     * @param value Array of values to be assigned
     * @param length Length of the value array
     * @param chain number of chain (starting from zero) to modify
     *
     * @see SArray#setValue
     */
    void setValue(double const *value, unsigned long length, unsigned int chain);
    /**
     * Indicates whether a node is discrete-valued or not.
     * @see SArray#isDiscreteValued
     */
    virtual bool isDiscreteValued() const = 0;
    /**
     * Returns a pointer to the start of the array of values for 
     * the given chain.
     */
    double const *value(unsigned int chain) const;
    /**
     * Returns the length of the value array
     */
    unsigned long length() const;
    /**
     * Returns the dimensions of the Node
     */
    std::vector<unsigned long> const &dim() const;
    /**
     * Swaps the values in the given chains
     */
    void swapValue(unsigned int chain1, unsigned int chain2);

    void addChild(StochasticNode *node) const;
    void removeChild(StochasticNode *node) const;
    void addChild(DeterministicNode *node) const;
    void removeChild(DeterministicNode *node) const;
    virtual void unlinkParents() = 0;
	
    /**
     * Returns the log of the density of a StochasticNode
     * given the current parameter values. For a ConstantNode
     * or a DeterministicNode this will simply return 0.
     *
     * @param chain Number of chain (starting from zero) for which
     * to evaluate log density.
     *
     * @param type Indicates whether the full probability density
     * function is required (PDF_FULL) or whether partial calculations
     * are permitted (PDF_PRIOR, PDF_LIKELIHOOD). See PDFType for
     * details.
     */
    virtual double logDensity(unsigned int chain, PDFType type) const = 0;
    virtual double KL(unsigned int chain1, unsigned int chain2, RNG *rng,
		      unsigned int nrep) const = 0;
	
    /**
     * Used by dumpNodeNames to gather a specific subset of node types:
     */
    virtual bool isConstant() const = 0;
    virtual bool isDeterministic() const = 0;
    virtual bool isStochastic() const = 0;
	
};

/**
 * Calculate the number of chains of parameter vector. Returns 0
 * if the parameters are inconsistent
 */
unsigned int countChains(std::vector<Node const *> const &parameters);

} /* namespace jags */

#endif /* NODE_H_ */
