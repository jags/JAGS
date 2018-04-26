#ifndef METROPOLIS_H_
#define METROPOLIS_H_

#include <sampler/MutableSampleMethod.h>
#include <vector>

namespace jags {

class StochasticNode;

/**
 * @short Metropolis-Hastings sampling method
 *
 * This class is used by Metropolis Hastings samplers.  It provides
 * only basic infrastructure.  
 *
 * The Metropolis class provides no update member function. A subclass
 * of Metropolis must provide this. It should contain one or more
 * calls to Metropolis#setValue, calculate the acceptance probability,
 * and then call the function Metropolis#accept.
 */
class Metropolis : public MutableSampleMethod
{
    std::vector<double> _last_value;
    bool _adapt;
    Metropolis(Metropolis const &);
    Metropolis &operator=(Metropolis const &);
public:
    Metropolis(std::vector<double> const &value);
    ~Metropolis();
    /**
     * Gets the current value array of the Metropolis object. 
     */
    virtual void getValue(std::vector<double> &value) const = 0;
    /**
     * Sets the value of the Metropolis object. 
     *
     * @param value Vector of values 
     */
    virtual void setValue(std::vector<double> const &value) = 0;
    /**
     * Accept current value with probabilty prob. If the current value
     * is not accepted, the Metropolis object reverts to the value at
     * the last successful call to accept.
     *
     * @param rng Random number generator.
     *  
     * @param prob Probability of accepting the current value.
     *
     * @returns success indicator
     */
    bool accept(RNG *rng, double prob);
    /**
     * Rescales the proposal distribution. This function is called by
     * Metropolis#accept when the sampler is in adaptive
     * mode. Rescaling may depend on the acceptance probability.
     *
     * @param prob Acceptance probability
     */
    virtual void rescale(double prob) = 0;
    /**
     * The Metropolis-Hastings method is adaptive. The process of
     * adaptation is specific to each subclass and is defined by the
     * rescale member function
     */
    bool isAdaptive() const;
    /**
     * Turns off adaptive mode
     */
    void adaptOff();
    /**
     * length of the value vector
     */
    unsigned long length() const;
};

} /* namespace jags */

#endif /* METROPOLIS_H_ */
