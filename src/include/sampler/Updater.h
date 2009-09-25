#ifndef UPDATER_H_ 
#define UPDATER_H_

#include <vector>
#include <string>
#include <set>

class StochasticNode;
class DeterministicNode;
class Node;
class Graph;
class RNG;

/**
 * @short Updates a set of stochastic nodes
 *
 * The Updater class is a helper class used by a Sampler. It allows
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
 * A vector of nodes in an acyclic Graph is in "forward sampling
 * order" if node B always appears after node A when there is a path
 * from A to B. Note that forward sampling order is not uniquely
 * determined.
 */
class Updater {
  unsigned int _length;
  std::vector<StochasticNode *> _nodes;
  std::vector<StochasticNode const *> _stoch_children;
  std::vector<DeterministicNode*> _determ_children;
public:
  /**
   * Constructs an Updater for the given vector of nodes.  
   *
   * @param nodes Vector of Nodes to be sampled 
   *
   * @param graph  Graph within which sampling is to take place. It is
   * an error if this Graph does not contain all of the Nodes to be sampled.
   */
  Updater(std::vector<StochasticNode *> const &nodes, Graph const &graph);
  /**
   * Constructs an Updater for a single node
   */
  Updater(StochasticNode * node, Graph const &graph);
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
  void setValue(double const * value, unsigned int length, unsigned int chain)
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
  std::vector<StochasticNode const*> const &stochasticChildren() const;
  /**
   * Returns the immediate deterministic descendendants of the sampled
   * nodes, in forward sampling order
   */
  std::vector<DeterministicNode*> const &deterministicChildren() const;
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
   * Static function that identifies the Marginal Stochastic Children
   * and the Immediate Deterministic Descendants of the given nodes
   * within the given graph.
   *
   * @param nodes Set of Nodes whose descendants are to be classified.
   *
   * @param graph Graph within which calculations are to take place.
   * Nodes outside of this graph will be ignored.
   *
   * @param stoch_nodes Vector which will contain the Marginal
   * Stochastic Children on exit.
   *
   * @param dtrm_nodes Vector which will contain the Immediate 
   * Deterministic Descendants, in forward sampling order, on exit.
   */
  static void classifyChildren(std::vector<StochasticNode *> const &nodes,
			       Graph const &graph,
			       std::vector<StochasticNode const*> &stoch_nodes,
			       std::vector<DeterministicNode*> &dtrm_nodes);
  /**
   * Simplified version of classifyChildren that finds only the stochastic
   * children within the graph.
   *
   * @param nodes Set of Nodes whose stochastic children are to be found.
   *
   * @param graph Graph within which calculations are to take place.
   * Nodes outside of this graph will be ignored.
   *
   * @param children Set which will contain the marginal stochastic
   * children on exit.
   */
  static void getStochasticChildren(std::vector<StochasticNode *> const &nodes,
				    Graph const &graph,
				    std::set<StochasticNode const*> &children);
};

unsigned int nchain(Updater const *updater);

#endif /* UPDATER_H_ */
