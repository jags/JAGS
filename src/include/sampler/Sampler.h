#ifndef SAMPLER_H_
#define SAMPLER_H_

#include <vector>
#include <string>
#include <set>

namespace jags {

struct RNG;
class StochasticNode;
class GraphView;

/**
 * @short Updates a set of stochastic nodes
 *
 * A sampler updates a set of stochastic nodes at each iteration.
 * In addition, the sampler is responsible for updating the immediate
 * deterministic descendants of the sampled nodes.
 *
 * The Sampler class is designed to be as abstract as possible to give
 * maximum freedom in implementing samplers.  In practice, all the
 * Samplers in the modules distributed with JAGS are implementations
 * of one of two subclasses -- MutableSampler or ImmutableSampler --
 * depending on whether the Sampler does or does not have a state that
 * it must preserve between iterations.
 */
class Sampler {
    GraphView *_gv;
public:
    /**
     * Constructor
     *
     * @param gv GraphView providing a view of the sampling graph for
     * the sampler. This should be dynamically allocated by the
     * SamplerFactory: the Sampler takes ownership of the GraphView
     * and deletes it when the destructor is called.
     */
    Sampler(GraphView *gv);
    virtual ~Sampler();
    /**
     * Returns the vector of stochastic nodes sampled by the Sampler
     */
    std::vector<StochasticNode*> const &nodes() const;
    /**
     * Updates the sampled nodes and their immediate deterministic
     * descendants for the given chain.
     *
     * @param chain Chain number
     * @param rng Pseudo-random number generator.
     */
    virtual void update(unsigned int chain, RNG * rng) = 0;
    /**
     * When a sampler is constructed, it may be in adaptive mode, which
     * allows it to adapt its behaviour for increased
     * efficiency. However, a sampler in adaptive mode may not converge
     * to the correct target distribution. This function turns off
     * adaptive mode, so that valid samples can be collected from the
     * sampler.
     */
    virtual void adaptOff() = 0;
    /**
     * The adaptOff function may be called at any time. Premature
     * ending of adaptive mode may result in an extremely inefficient
     * sampler.  Therefore the checkAdaptation function implements an
     * efficiency test that returns true if it is safe to call the
     * adaptOff function.  Samplers that have no adaptive mode should
     * simply return true.
     */
    virtual bool checkAdaptation() const = 0;
    /**
     * Indicates whether the sampler has an adaptive mode.
     */
    virtual bool isAdaptive() const = 0;
    /**
     * Returns a name for the sampler which should describe the method
     * it uses to update the nodes.
     */
    virtual std::string name() const = 0;
};

} /* namespace jags */

#endif /* SAMPLER_H_ */
