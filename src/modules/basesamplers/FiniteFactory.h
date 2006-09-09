#ifndef FINITE_FACTORY_H_
#define FINITE_FACTORY_H_

#include <sampler/GibbsFactory.h>
class StochasticNode;
class Graph;

/**
 * @short Factory object for finite samplers
 */
class FiniteFactory : public GibbsFactory
{
 public:
  bool canSample(StochasticNode *snode, Graph const &graph) const;
  Sampler *makeGibbsSampler(StochasticNode *snode, Graph const &graph,
                            unsigned int nchain) const;
};

#endif /* FINITE_FACTORY_H_ */
