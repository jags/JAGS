#ifndef SAMPLER_H_
#define SAMPLER_H_

#include <vector>
#include <set>

class StochasticNode;
class DeterministicNode;
class Node;
class Graph;
class RNG;

/**
 * A sampler updates a set of stochastic nodes.  It is also
 * responsible for updating the immediate deterministic descendants of
 * those nodes (see below for a definition).  Sampling takes place in
 * the context of a Graph, which must contain the sampled nodes. Any
 * descendents of these nodes outside of the Graph are ignored when
 * updating.
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
 * A vector of nodes is in "forward sampling order" if node B always
 * appears after node A when there is a path from A to B. Note that
 * forward sampling order is not uniquely determined.
 */
class Sampler {
protected:
  const unsigned int _chain;
private:
  std::vector<StochasticNode *> _nodes;
  std::vector<StochasticNode const *> _stoch_children;
  std::vector<DeterministicNode*> _determ_children;
public:
  /**
   * Constructs a sampler for the given vector of nodes.  
   *
   * @param nodes Vector of Nodes to be sampled 
   * @param graph  Graph within which sampling is to take place. It is
   * an error if this Graph does not contain all of the Nodes to be sampled.
   * @param chain  Index of chain sampled by this sampler.
   */
  Sampler(std::vector<StochasticNode *> const &nodes, Graph const &graph,
	  unsigned int chain);
  virtual ~Sampler();
  /**
   * Returns the vector of sampled nodes
   */
  std::vector<StochasticNode *> const &nodes() const;
  /**
   * Sets the values of the sampled nodes.  The immediate
   * deterministic descendants are automatically updated.
   * This function should be called by the update function.
   *
   * @param value vector of pointers to the beginning of an array of
   * values for each node.  
   * @param length vector of lengths for the corresponding arrays.
   */
  void setValue(std::vector<double const *> const &value, 
		std::vector<unsigned long> const &length);
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
   * Returns the chain number that is sampled by this sampler
  */
  unsigned int chain();
  /**
   * Calculates the log conditional density of the sampled nodes,
   * given all other nodes in the graph that was supplied to the
   * constructor, plus the parents of the nodes (which may be outside
   * the graph).  The log full conditional is calculated up to an
   * additive constant.
   */
  double logFullConditional();
  /**
   * Every sampler must update the vector of nodes and its immediate
   * deterministic descendants using the update function.
   *
   * @param rng Pseudo-random number generator function.
   */
  virtual void update(RNG *rng) = 0;
  /**
   * When a sampler is constructed, it may be in "burnin" mode, which
   * allows it to adapt its behaviour for increased efficiency. However,
   * a sampler in burnin mode may not necessarily supply valid updates.
   *
   * The purpose of this function is to turn off burnin mode, so that
   * valid samples can be collected from the sampler.
   */
  virtual void burninOff() = 0;
  /**
   * Static function that identifies the Marginal Stochastic Children
   * and the Immediate Deterministic Descendants of the given nodes
   * within the given graph.
   *
   * @param nodes Set of Nodes whose descendants are to be classified.
   * @param graph Graph within which calculations are to take place.
   * Nodes outside of this graph will be ignored.
   * @param stoch_nodes Empty vector which will contain the Marginal
   * Stochastic Children on exit.
   * @param dtrm_nodes Empty vector which will contain the Immediate 
   * Deterministic Descendants, in forward sampling order, on exit.
   */
  static void classifyChildren(std::vector<StochasticNode *> const &nodes,
			       Graph const &graph,
			       std::vector<StochasticNode const*> &stoch_nodes,
			       std::vector<DeterministicNode*> &dtrm_nodes);
};

#endif /* SAMPLER_H_ */
