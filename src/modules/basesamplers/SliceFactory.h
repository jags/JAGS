#ifndef SLICE_FACTORY_H_
#define SLICE_FACTORY_H_

#include <sampler/GibbsFactory.h>
class StochasticNode;
class Graph;

/**
 * @short Factory object for slice samplers
 */
class SliceFactory : public GibbsFactory
{
 public:
  bool canSample(StochasticNode *snode, Graph const &graph) const;
  Sampler *makeGibbsSampler(StochasticNode *snode, Graph const &graph, unsigned int chain) const;
};

#endif /* SLICE_FACTORY_H_ */
