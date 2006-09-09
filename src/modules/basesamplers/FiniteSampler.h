#ifndef FINITE_SAMPLER_H_
#define FINITE_SAMPLER_H_

#include <sampler/GibbsSampler.h>

class StochasticNode;

/**
 * Sampler for discrete distributions with support on a finite set.
 */
class FiniteSampler : public GibbsSampler {
public:
    FiniteSampler(StochasticNode *node, Graph const &graph, unsigned int chain);
    ~FiniteSampler();

    void update(RNG *rng);
    static bool canSample(StochasticNode const *node, Graph const &graph);
    void burninOff();
};

#endif /* FINITE_SAMPLER_H_ */
