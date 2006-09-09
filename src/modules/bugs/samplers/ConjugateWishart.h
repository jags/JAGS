#ifndef CONJUGATE_WISHART_H_
#define CONJUGATE_WISHART_H_

#include "ConjugateSampler.h"

/**
 * Conjugate sampler for the Wishart distribution
 *
 * A Conjugate Wishart sampler for stochastic node "snode" can be
 * created if snode has a Wishart distribution and its 
 * stochastic children are multivariate normal with precision matrix
 * snode. The mean of the children may not depend on snode.
 */
class ConjugateWishart : public ConjugateSampler {
  bool _fast;
public:
  ConjugateWishart(StochasticNode *snode, Graph const &graph, unsigned int chain);
  ~ConjugateWishart();
  
  static bool canSample(StochasticNode *snode, Graph const &graph);
  void update(RNG *);
};

#endif /* CONJUGATE_WISHART_H_ */
