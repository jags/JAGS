#ifndef CONJUGATE_DIRICHLET_H_
#define CONJUGATE_DIRICHLET_H_

#include "ConjugateSampler.h"

/**
 * @short Conjugate sampler for Dirichlet distribution.  
 *
 * A conjugate sampler may be created for the stochastic node "snode"
 * if it has a Dirichlet distribution and its marginal stochastic
 * children have either a categorical distribution (dcat) or a
 * multinomial distribution (dmult) with snode as the probability
 * vector.  In the case of a multinomial distribution, the sample size
 * must not depend on snode.  The immediate deterministic children
 * must be either identity functions or mixture nodes.
 */
class ConjugateDirichlet : public ConjugateSampler {
 public:
  ConjugateDirichlet(StochasticNode *snode, Graph const &graph, unsigned int chain);
  ~ConjugateDirichlet();
  
  static bool canSample(StochasticNode *snode, Graph const &graph);
  void update(RNG *rng);
};

#endif /* CONJUGATE_DIRICHLET_H_ */
