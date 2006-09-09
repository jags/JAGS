#ifndef CENSORED_SAMPLER_H_
#define CENSORED_SAMPLER_H_

#include <sampler/GibbsSampler.h>

/**
 *
 * A sampler for censored real-valued distributions.
 *
 */
class Censored : public GibbsSampler {
 public:
  Censored(StochasticNode *snode, Graph const &graph, unsigned int chain);
  ~Censored();

  static bool canSample(StochasticNode *snode, Graph const &graph);
  void update(double (*unif_rand)(void));
  void burninOff();
};

#endif /* CONJUGATE_NORMAL_H_ */
