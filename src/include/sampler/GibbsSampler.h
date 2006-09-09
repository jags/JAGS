#ifndef GIBBS_SAMPLER_H_
#define GIBBS_SAMPLER_H_

#include <sampler/Sampler.h>

#include <set>

class StochasticNode;
class Node;

/**
 * A Gibbs sampler acts on a single stochastic node (the node may
 * itself be multivariate). The GibbsSampler class is abstract because
 * it does not implement Sampler##update, but provides convenience
 * functions to make this easier when there is only one node in the
 * sampler.
 */
class GibbsSampler : public Sampler {
 public:
  /** Construct a sampler for the given node.
   * @see Sampler::Sampler
   */
  GibbsSampler(StochasticNode *node, Graph const &graph, unsigned int chain);
  /** Returns the node being sampled */
  StochasticNode *node();
  /** Sets the value of the node and updates its immediate
   * deterministic descendants
   * @see Sampler::setValue
   */
  void setValue(double const *value, unsigned long length);
};

#endif /* GIBBS_SAMPLER_H_ */
