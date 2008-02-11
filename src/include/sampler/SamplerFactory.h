#ifndef SAMPLER_FACTORY_H_
#define SAMPLER_FACTORY_H_

#include <vector>
#include <set>

class Sampler;
class StochasticNode;
class Graph;

/**
 * @short Factory for Sampler objects
 */
class SamplerFactory
{
public:
    virtual ~SamplerFactory();
    /**
     * Finds nodes in the set of stochastic nodes that can be sampled
     * within the given graph, removes them from the set and adds a
     * sampler to the vector of samplers.
     */
    virtual void 
	makeSampler(std::set<StochasticNode*> &nodes, Graph const &graph,
		    std::vector<Sampler*> &samplers) const = 0;
};

#endif /* SAMPLER_FACTORY_H_ */
