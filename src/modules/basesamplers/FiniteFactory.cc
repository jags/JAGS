#include <config.h>
#include "FiniteSampler.h"
//#include <graph/StochasticNode.h>
#include "FiniteFactory.h"

bool FiniteFactory::canSample(StochasticNode * snode, Graph const &graph) const
{
  return FiniteSampler::canSample(snode, graph);
}

Sampler *FiniteFactory::makeGibbsSampler(StochasticNode *snode,
					 Graph const &graph,
                                         unsigned int chain) const
{
  return new FiniteSampler(snode, graph, chain);
}
