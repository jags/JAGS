#ifndef MODEL_H_
#define MODEL_H_

#include <graph/Graph.h>

#include <vector>
#include <list>

class Sampler;
class TraceMonitor;
class SamplerFactory;
class RNG;
class RNGFactory;

struct ChainInfo
{
  RNG *rng;
  std::vector<Sampler*> samplers;
  std::list<TraceMonitor*> monitors;
  unsigned long iteration;
};

/**
 * @short Graphical model 
 *
 * The purpose of the model class is to collect together all the
 * elements necessary to run a MCMC sampler on a graphical model.
 */
class Model {
  unsigned int _nchain;
  std::vector<ChainInfo> _chain_info;
  Graph _graph;
  std::vector<Node*> _nodes;
  std::set<Node*> _extra_nodes;
  std::vector<Node*> _sampled_extra;
  std::set<Node*>  _monitored_nodes;
  bool _is_graph_checked;
  bool _is_initialized;
  bool _can_sample;
  void chooseRNGs();
public:
  Model(unsigned int nchain);
  virtual ~Model();

  Graph &graph();
  /**
   * Checks that the graph is closed and acyclic. A runtime_error
   * is thrown if it is not
   */
  void checkGraph();
  /**
   * Returns true if the graph has been checked
   */
  bool isGraphChecked();
  /**
   * Initializes the model.  
   *
   * Firstly, the Node##initialize function is called for all nodes in
   * forward sampling order.
   * 
   * Secondly, any chain without a random number generator is assigned
   * an RNG object by traversing the list of RNGFactories.
   *
   * @see Node##initialize, Model##rngFactories
   */
  void initialize();
  /** Returns true if the model has been initialized */
  bool isInitialized();
  /*
   * Selects samplers. For each chain, samplers are selected by
   * traversing the list of SamplerFactories in order. If there are
   * any informative stochastic nodes left without samplers after all
   * factories have been tried, then a runtime error is thrown
   *
   * @see Model##samplerFactories
   */
  void chooseSamplers();
  /**
   * Returns true if chooseSamplers has been called.
   */
  bool canSample() const;
  /** 
   * Updates the model by the given number of iterations. A
   * logic_error is thrown if the model is uninitialized.
   *
   * @param niter Number of iterations to run
   */
  void update(long niter);
  /**
   * Returns the current iteration number for the given chain. Chains
   * are enumerated starting from 0.
   */
  long iteration(unsigned int chain) const;
  /**
   * Creates a monitor for the given node, with given thinning
   * interval If the node is already being monitored, a NodeError is
   * thrown.  A side effect is to turn off burnin mode for all
   * samplers in the model.
   */
  void setMonitor(Node *node, int thin);
  /**
   * Clears the monitor associated with the given node. If there is no
   * such monitor, then the function has no effect.
   */
  void clearMonitor(Node *node);
  /**
   * Returns the list of monitors
   */
  std::list<TraceMonitor*> const &monitors(unsigned int chain) const;
  /**
   * After the model is initialized, extra uninformative nodes may be
   * added to the graph (For example, a node may be added that
   * calculates the deviance of the model).  The model takes
   * responsibility for updating the extra node.
   *
   * The extra node cannot be observed, it must not already be in the
   * model graph, it may not have any children, and all of its parents
   * must be in the graph.
   */
  void addExtraNode(Node *node);
  /**
   * Access the list of sampler factories, which is common to all
   * models. This is used during initialization to choose samplers.
   *
   * @seealso Model##chooseSamplers
   */
  static std::list<SamplerFactory const *> &samplerFactories();
  /**
   * Access the list of RNG factories, which is common to all models
   */
  static std::list<RNGFactory *> &rngFactories();
  /**
   * Returns the number of chains in the model
   */
  unsigned int nchain() const;
  /**
   * Returns the RNG object associated with the given chain. If no RNG
   * has been assigned, this a NULL pointer is returned.
   */
  RNG *rng(unsigned int nchain) const;
  /**
   * Assigns a new RNG object to the given chain. The list of
   * RNGFactory objects is traversed and each factor is requested to
   * generate a new RNG object of the given name.
   *
   * @return success indicator.
   */
  bool setRNG(std::string const &name, unsigned int chain);
  /**
   * Assigns an existing RNG object to the given chain
   *
   * @return success indicator
   */
  bool setRNG(RNG *rng, unsigned int chain);
};

#endif /* MODEL_H_ */
