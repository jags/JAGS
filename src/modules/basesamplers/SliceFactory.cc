#include <config.h>

#include "RealSliceSampler.h"
#include "DiscreteSliceSampler.h"
#include "SliceFactory.h"

#include <stdexcept>

using std::logic_error;

bool SliceFactory::canSample(StochasticNode * snode, Graph const &graph) const
{
    return RealSliceSampler::canSample(snode, graph) || 
	DiscreteSliceSampler::canSample(snode, graph);
}

Sampler *SliceFactory::makeGibbsSampler(StochasticNode *snode,
					Graph const &graph,
                                        unsigned int chain) const
{
    if (RealSliceSampler::canSample(snode, graph)) {
	return new RealSliceSampler(snode, graph, chain);
    }
    else if (DiscreteSliceSampler::canSample(snode, graph)) {
	return new DiscreteSliceSampler(snode, graph, chain);
    }
    else {
	throw logic_error("Unable to construct slice sampler");
    }
}
