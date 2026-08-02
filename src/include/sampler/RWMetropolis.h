#ifndef RW_METROPOLIS_H_
#define RW_METROPOLIS_H_

#include <sampler/Metropolis.h>
#include <sampler/StepAdapter.h>

namespace jags {

/**
 * @short Random Walk Metropolis-Hastings sampling method
 *
 * This class provides an update function which modifies the current
 * value by a random walk. This step size is adapted to achieve the
 * target acceptance rate using a noisy gradient algorithm.
 */
class RWMetropolis : public Metropolis
{
    StepAdapter _step_adapter;
    double _pmean;
    unsigned int _niter;
public:
    /**
     * Constructor. 
     *
     * @param value Initial value vector.
     * @param step Initial step size for the random walk updates.
     * @param accept See StepAdapter
     * @param delta See StepAdapter
     * @param nstart See StepAdapter
     * @param min_step See StepAdapter
     */
    RWMetropolis(std::vector<double> const &value, double step,
		 double accept = 0.234, double delta = 2.17, double nstart = 28, double min_step = 0);
    ~RWMetropolis() override;
    /**
     * Updates the current value by adding a random increment.
     */
    void update(RNG *rng) override;
    /**
     * Modifies the step size to achieve the target acceptance
     * probability using a noisy gradient algorithm
     *
     * @param prob acceptance probability at current update
     */
    void rescale(double prob) override;
    /**
     * The RWMetropolis method keeps a running mean of the acceptance
     * probabilities, which is updated every time the rescale function
     * is called. The checkAdaptation function returns true if the logit
     * of the running mean is within 0.50 of the target.
     */
    bool checkAdaptation() const override;
    /**
     * Modifies the given value vector in place by adding an
     * independent normal increment to each element.  It can be
     * overridden to provide non-normal increments, or a random walk
     * on some transformed scale (but see RMetropolis#logJacobian).
     *
     * Note that this function does not modify the value of the
     * RWMetropolis object.
     * 
     */
    virtual void step(std::vector<double> &value, double s, RNG *rng) const;
    /**
     * If the random walk takes place on a transformed scale
     * (e.g. log, logistic), then the log density of the target
     * distribution must be penalized by the log Jacobian of the
     * transformation.
     *
     * This function calculates the log Jacobian at the given value.
     * By default, the random walk takes place on the original scale
     * and therefore the penalty is zero.
     */
    virtual double logJacobian(std::vector<double> const &value) const;
    /**
     * Returns the log of the probability density function of the target
     * distribution at the current value.
     */
    virtual double logDensity() const = 0;
};

} /* namespace jags */

#endif /* RW_METROPOLIS_H_ */
