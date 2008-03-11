#ifndef STOCHASTIC_NODE_H_
#define STOCHASTIC_NODE_H_

#include <graph/Node.h>

class Distribution;
class RNG;

/**
 * @short Node defined by the BUGS-language operator ~
 *
 * Stochastic nodes represent the random variables that are the
 * fundamental building blocks of a Bayesian hierarchical model. In
 * the BUGS language, they are defined on the left hand side of a
 * stochastic relation. For example, the relation 
 *
 * <pre>y ~ dnorm(mu, tau) T(L, U)</pre> 
 *
 * defines y to be a normally distributed random variable with parameters
 * mu,  tau, L, and U (mean, precision, lower bound, upper bound). The
 * last two parameters, defined by the T(,) construct, are optional. If
 * they are supplied, then the distribution of the node is truncated
 * to lie in the range (L, U). Not all distributions can be truncated.
 *
 * JAGS allows you to define stochastic nodes that are, in fact,
 * not random at all, but are deterministic functions of their parameters.
 * A common example is the dinterval distribution
 *
 * <pre>group[i] ~ dinterval(true[i], cutpoints[1:N])</pre>
 *
 * where the value of group[i] is determined by where the value of
 * true[i] falls in the vector of supplied cutpoints.  In this case,
 * the stochastic node leads a double life. If it is observed, then it
 * is considered a random variable, and generates a likelihood for its
 * stochastic parents.  If it is unobserved then it is treated as a
 * deterministic function of its parents, just as if it were a
 * LogicalNode.
 * 
 * @see Distribution
 */
class StochasticNode : public Node {
    Distribution const * const _dist;
    std::vector<std::vector<double const *> > _parameters;
    std::vector<std::vector<unsigned int> > _dims;
    Node const *_lower;
    Node const *_upper;
public:
    /**
     * Constructs a new StochasticNode given a distribution, a vector
     * of parent nodes, considered as parameteres to the distribution,
     * and, optionally, upper and lower bounds. If bounds are given
     * then the distribution of the constructed StochasticNode is
     * truncated at the value of the bounds. 
     */
    StochasticNode(Distribution const *dist, 
                   std::vector<Node const *> const &parameters,
                   Node const *lower=0, Node const *upper=0);
    ~StochasticNode();
    /**
     * Returns a pointer to the Node that defines the lower bound, if
     * the distribution is truncated, or a NULL pointer otherwise.
     */
    Node const *lowerBound() const;
    /**
     * Returns a pointer to the Node that defines the upper bound, if
     * the distribution is truncated, or a NULL pointer otherwise.
     */
    Node const *upperBound() const;
    /**
     * Returns a pointer to the value of the lower bound for the given
     * chain, if the distribution is truncated, or a NULL pointer
     * otherwise.
     */
    double const *lowerLimit(unsigned int chain) const;
    /**
     * Returns a pointer to the value of the upper bound for the given
     * chain, if the distribution is truncated, or a NULL pointer
     * otherwise.
     */
    double const *upperLimit(unsigned int chain) const;
    /**
     * Returns a pointer to the Distribution of the StochasticNode.
     */
    Distribution const *distribution() const;
    /**
     * Returns a vector of parameter values for the Distribution of
     * the stochastic node. Each element of the vector is a pointer
     * to the start of an array of doubles. It is assumed that these
     * arrays are of the correct size.
     *
     * @param chain Index number of the chain for which parameters are
     * requested.  
     */
    std::vector<double const *> const &parameters(unsigned int chain) const;
    /**
     * Returns a vector of dimensions for the parameters for the
     * distribution. These are the parameters of the parent Nodes
     * supplied to the constructor.
     */
    std::vector<std::vector<unsigned int> > const &parameterDims() const;
    /**
     * Returns the log of the prior density of the StochasticNode
     * given the current parameter values.
     */
    double logDensity(unsigned int chain) const;
    /**
     * Draws a random sample from the prior distribution of the node
     * given the current values of it's parents, and sets the Node
     * to that value.
     *
     * @param rng Random Number Generator object
     *
     * @param chain Index umber of chain to modify
     */
    void randomSample(RNG *rng, unsigned int chain);
    /**
     * A deterministic sample for a stochastic node sets it to a
     * "typical" value of the prior distribution, given the current
     * values of its parents. The exact behaviour depends on the
     * Distribution used to define the StochasticNode, but it will
     * usually be the prior mean, median, or mode.
     */
    void deterministicSample(unsigned int chain);
    /**
     * Ensures that the values of the stochastic node parents are valid
     * 
     * @see Distribution#checkParameterValue
     */
    bool checkParentValues(unsigned int chain) const;
    /**
     * Stochastic nodes are normally considered to represent random
     * variables in the model. Hence, this function usually returns
     * true. However, there is an important exception.
     * 
     * If the number of degrees of freedom of the node's Distribution
     * is zero, then the value of the node is a deterministic function
     * of it's parents. In this case, if the node is unobserved, it is 
     * not considered to be a random variable.
     */
    bool isRandomVariable() const;
    /**
     * Checks that the parameters are within the valid range determined
     * by the Distribution.
     *
     * @param chain Index number of chain to check
     *
     * @see Distribution#checkParameterValue
     */
    bool checkParameterValue(unsigned int chain) const;    
    /**
     * Stochastic nodes are never linear functions of their
     * parameters. This function always returns false.
     */
    bool isLinear(GraphMarks const &linear_marks, bool fixed) const;
    /**
     * Stochastic nodes are never scale functions. This function
     * always returns false.
     */
    bool isScale(GraphMarks const &scale_marks, bool fixed) const;
    std::string deparse(std::vector<std::string> const &parameters) const;
};

/**
 * Wrapper function that dynamically casts a Node pointer to a
 * StochasticNode pointer.
 */
StochasticNode const *asStochastic(Node const *node);

/**
 * Number of degrees of freedom of a node
 *
 * @see Distribution#df
 */
unsigned int df(StochasticNode const *snode);

/**
 * Indicates whether the distribution of the node is bounded
 * either above or below.
 */
bool isBounded(StochasticNode const *node);

/**
 * Writes the lower and upper limits of the support of a given
 * stochastic node to the supplied arrays. If the node has upper and
 * lower bounds then their values are taken into account in the
 * calculation.
 *
 * @param lower pointer to start of an array that will hold the lower 
 * limit of the support
 *
 * @param lower pointer to start of an array that will hold the upper 
 * limit of the support
 *
 * @param length size of the lower and upper arrays.
 *
 * @param node Stochastic node to query
 *
 * @param chain Index number of chain to query
 *
 * @see Distribution#support
 */
void support(double *lower, double *upper, unsigned int length,
             StochasticNode const *node, unsigned int chain);

/**
 * Returns true if the upper and lower limits of the support of the
 * stochastic node are fixed. Upper and lower bounds are taken into account
 *
 * @see Distributin#isSupportFixed
 */
bool isSupportFixed(StochasticNode const *node);

#endif /* STOCHASTIC_NODE_H_ */

