#ifndef RW_DSUM_H_
#define RW_DSUM_H_

#include <sampler/Metropolis.h>
#include <sampler/StepAdapter.h>

class Graph;
class GraphView;

namespace jags {
namespace bugs {

/**
 * @short Random Walk Metropolis-Hastings for parents of DSum nodes
 */
class RWDSum : public Metropolis
{
    GraphView const *_gv;
    unsigned int _chain;
    StepAdapter _step_adapter;
    double _pmean;
    unsigned int _niter;
protected:
    StochasticNode const *_dsnode;
public:
    /**
     * Constructor. 
     *
     * @param step Initial step size for the random walk updates.
     */
    RWDSum(GraphView const *gv, unsigned int chain, double step);
    /**
     * Updates the current value by adding a random increment.
     */
    void update(RNG *rng);
    /**
     * Modifies the step size to achieve the target acceptance
     * probability using a noisy gradient algorithm
     *
     * @param prob acceptance probability at current update
     */
    void rescale(double prob);
    /**
     * The RWDSum method keeps a running mean of the acceptance
     * probabilities, which is updated every time the rescale function
     * is called. The checkAdaptation function returns true if the logit
     * of the running mean is within 0.50 of the target.
     */
    bool checkAdaptation() const;
    /**
     * Does a random walk step.  Note that this function does not
     * modify the value of the RWDSum object.
     *
     * @param value Value vector. It represents an nrow by ncol matrix
     *        of values stored in column order.
     * @param nrow Number of rows of the value matrix
     * @param ncol Number of columns of the value matrix
     * @param s Step size
     * @param rng Random number generator
     */
    virtual void step(std::vector<double> &value, unsigned long nrow,
		      unsigned long ncol, double s, RNG *rng) const = 0;
    /**
     * Tests whether nodes can be sampled
     *
     * @param nodes Vector of stochastic nodes to be sampled
     * @param graph Graph within which sampling takes place
     * @param discrete Flag indicating whether the nodes are discrete-valued
     *        or not
     */
    static bool canSample(std::vector<StochasticNode *> const &nodes,
			  Graph const &graph, bool discrete, bool multinom);
    void setValue(std::vector<double> const &value);
    void getValue(std::vector<double> &value) const;

};

}}

#endif /* RW_DSUM_H_ */
