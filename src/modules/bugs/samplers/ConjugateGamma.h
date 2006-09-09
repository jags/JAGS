#ifndef CONJUGATE_GAMMA_H_
#define CONJUGATE_GAMMA_H_

#include "ConjugateSampler.h"

class ConjugateGamma : public ConjugateSampler {
    double *_coef;
    void calCoef();
 public:
  ConjugateGamma(StochasticNode *snode, Graph const &graph, unsigned int chain);
  ~ConjugateGamma();
	
  static bool canSample(StochasticNode *snode, Graph const &graph);
  void update(RNG *rng);
};

#endif /* CONJUGATE_GAMMA_H_ */
