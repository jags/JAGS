#ifndef CONJUGATE_NORMAL_H_
#define CONJUGATE_NORMAL_H_

#include "ConjugateSampler.h"

/**
 * Conjugate sampler for the normal distribution.
 * A conjugate sampler can be created for a stochastic node "snode"
 * if it has a normal distribution, all of its marginal stochastic
 * children are normal, they depend only on snode through the mean
 * (and not the precision) and the mean is a linear function of
 * snode.
 */
class ConjugateNormal : public ConjugateSampler {
    double *_betas;
    unsigned long _length_betas;
    void calBeta();
 public:
  ConjugateNormal(StochasticNode *snode, Graph const &graph, unsigned int chain);
  ~ConjugateNormal();

  static bool canSample(StochasticNode *snode, Graph const &graph);
  void update(RNG *rng);
};

#endif /* CONJUGATE_NORMAL_H_ */
