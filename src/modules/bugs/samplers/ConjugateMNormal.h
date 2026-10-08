#ifndef CONJUGATE_M_NORMAL_H_
#define CONJUGATE_M_NORMAL_H_

#include "ConjugateMethod.h"

namespace jags {
    
    class Graph;

namespace bugs {

    class LMDesign;

/**
 * Conjugate sampler for multivariate normal distribution
 *
 * ConjugateMNormal is currently rather limited in scope.  It can
 * sample a multivariate normal node "snode" if all the stochastic
 * children of snode in the graph are themselves multivariate normal,
 * with mean snode, and where the precision matrix does not depend on
 * snode.
 *
 * This can (and should) be extended to allow marginal stochastic
 * children that are normally distribution and where the mean of the
 * children is a linear function of (some element of) snode.
 */
class ConjugateMNormal : public ConjugateMethod {
  double *_betas;
  LMDesign *_design;  
  const unsigned int _length_betas;
 public:
  ConjugateMNormal(SingletonGraphView const *gv);
  ~ConjugateMNormal() override;
  void update(unsigned int chain, RNG *rng) const override;
  static bool canSample(StochasticNode *snode, Graph const &graph);
};

}}

#endif /* CONJUGATE_M_NORMAL_H_ */
