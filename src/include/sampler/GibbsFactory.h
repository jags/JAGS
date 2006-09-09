#ifndef GIBBS_FACTORY_H_
#define GIBBS_FACTORY_H_

#include <sampler/SamplerFactory.h>
/**
 * @short Factory object for Samplers that sample a single node
 */
class GibbsFactory : public SamplerFactory
{
  /**
   * Determines whether the factory can produce a GibbsSampler for
   * the given node, within the given graph.
   */
  virtual bool canSample(StochasticNode *node, Graph const &graph) 
      const = 0;
  /**
   * Returns a dynamically allocated Sampler for a given node.
   */
  virtual Sampler *makeGibbsSampler(StochasticNode *node,
				    Graph const &graph,
                                    unsigned int chain) const = 0;
  /**
   * Implementation of SamplerFactory##makeSampler. This traverses
   * the graph, creating GibbsSampler, when possible, for each
   * individual node
   */
  void makeSampler(std::set<StochasticNode*> &nodes, 
                   Graph const &graph,
		   std::vector<std::vector<Sampler*> > &samplers) const;
};

#endif /* GIBBS_FACTORY_H */
