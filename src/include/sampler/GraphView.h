#ifndef GRAPH_VIEW_H_ 
#define GRAPH_VIEW_H_

#include <vector>
#include <string>
#include <set>

namespace jags {

class StochasticNode;
class DeterministicNode;
class Node;
class Graph;
struct RNG;

/**
 * @short Interface to a graph for Samplers.
 *
 * The GraphView class is a helper class used by a Sampler. It allows
 *  new values to be assigned to a vector of stochastic nodes.  When
 *  new values are assigned, it also updates the immediate
 *  deterministic descendants of those nodes (see below for a
 *  definition).  Sampling takes place in the context of a Graph,
 *  which must contain the sampled nodes. Any descendents of these
 *  nodes outside of the Graph are ignored when updating.
 *
 * Some terminology:
 *
 * The "immediate deterministic descendants" of a set of stochastic
 * nodes S are the descendants of S in the graph where all stochastic
 * nodes except those in S have been removed.
 *
 * The "marginal stochastic children" of a set of stochastic nodes S
 * are the children of S in the graph where all deterministic nodes
 * have been marginalized out.
 *
 * A vector of nodes in an acyclic Graph is in topological order if
 * node B always appears after node A when there is a path from A to
 * B. Note that topological order is not uniquely determined.
 */
class GraphView {
  unsigned int _length;
  std::vector<StochasticNode *> _nodes;
  std::vector<StochasticNode *> _stoch_children;
  std::vector<DeterministicNode*> _determ_children;
  bool _multilevel;
  void classifyChildren(std::vector<StochasticNode *> const &nodes,
			Graph const &graph,
			std::vector<StochasticNode *> &stoch_nodes,
			std::vector<DeterministicNode*> &dtrm_nodes, 
			bool allow_multilevel);
public:
  /**
   * Constructs a GraphView for the given vector of nodes.  
   *
   * @param nodes Vector of Nodes to be sampled 
   *
   * @param graph Graph within which sampling is to take place. It is
   *        an error if this Graph does not contain all of the Nodes
   *        to be sampled.
   *
   * @param allow_multilevel Indicates whether the GraphView may be
   *        multilevel, i.e.  some sampled nodes are stochastic
   *        children of other sampled nodes.  Many update methods do
   *        not work with multi-level models so this parameter is false
   *        by default.  If an attempt is made to create a multi-level graph
   *        view when the multilevel is false, a logic_error is thrown.
   */
  GraphView(std::vector<StochasticNode *> const &nodes, Graph const &graph,
	    bool allow_multilevel=false);
  /**
   * Returns the vector of sampled nodes.
   */
  std::vector<StochasticNode *> const &nodes() const;
  /**
   * Sets the values of the sampled nodes.  Their immediate
   * deterministic descendants are automatically updated.
   *
   * @param value Array of concatenated values to be applied to the 
   * sampled nodes.
   *
   * @param length Length of the value array. This must be equal to the
   * sum of the  lengths of the sampled nodes.
   *
   * @param chain Number of the chain (starting from zero) to be modified.
   */
  void setValue(double const * value, unsigned long length, unsigned int chain)
      const;
  void setValue(std::vector<double> const &value, unsigned int chain) const;
  void getValue(std::vector<double> &value, unsigned int chain) const;
  /**
   * Returns the total length of the sampled nodes.
   */
  unsigned int length() const;
  /**
   * Returns the marginal stochastic children of the sampled nodes.
   */
  std::vector<StochasticNode *> const &stochasticChildren() const;
  /**
   * Returns the immediate deterministic descendendants of the sampled
   * nodes, in topological order
   */
  std::vector<DeterministicNode*> const &deterministicChildren() const;
  /**
   * Tests whether the node depends deterministically on any of the
   * sampled nodes.  This function may be used by SamplerFactory
   * objects to test the validity of sampling methods. For example, we
   * may need to check that the mean, but not the variance of a
   * stochastic child depends on the sampled nodes.
   *
   * @param node Node to test
   *
   * @return true if node is in either the vector of sampled nodes or the
   * vector of deterministic children, false otherwise.
   */
  bool isDependent(Node const *node) const;
  /**
   * Calculates the log conditional density of the sampled nodes,
   * given all other nodes in the graph that was supplied to the
   * constructor, plus the parents of the nodes (which may be outside
   * the graph).  The log full conditional is calculated up to an
   * additive constant.
   *
   * @param chain Number of the chain (starting from zero) to query.
   */
  double logFullConditional(unsigned int chain) const;
  /**
   * Calculates the log prior density of the sampled nodes, i.e. the
   * density conditioned only on the parents.
   */
  double logPrior(unsigned int chain) const;
  /**
   * Calculates the log likelihood, which is added to the log prior
   * to give the log full conditional density
   */
  double logLikelihood(unsigned int chain) const;
  /**
   * Checks that the log density is finite for all sampled nodes and
   * all stochastic children. If any log density value is negative
   * infinite a NodeError exception is thrown.
   *
   * This utility function should be called by the constructor of any
   * sample method that needs to ensure the graph is in a valid state
   * before sampling.
   */
  void checkFinite(unsigned int chain) const;
};

unsigned int nchain(GraphView const *gv);

} /* namespace jags */

#endif /* GRAPH_VIEW_H_ */
