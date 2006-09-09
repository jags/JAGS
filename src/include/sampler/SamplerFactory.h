#ifndef SAMPLER_FACTORY_H_
#define SAMPLER_FACTORY_H_

#include <vector>
#include <set>

class Sampler;
class StochasticNode;
class Graph;

/**
 * @short Factory object for samplers
 */
class SamplerFactory
{
 public:
  virtual ~SamplerFactory() {};
  /**
   * Finds nodes in the set of stochastic nodes that can be sampled,
   * removes them from the set and adds one or more samplers to the
   * vector of samplers.
   */
  virtual void makeSampler(std::set<StochasticNode*> &nodes, 
			   Graph const &graph,
			   std::vector<std::vector<Sampler*> > &samplers)
  const = 0;
};

#endif /* SAMPLER_FACTORY_H_ */
