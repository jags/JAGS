#ifndef RNG_H_
#define RNG_H_

#include <string>
#include <vector>

namespace jags {

/**
 * Abstract class for a psuedo random number generator
 *
 * @short Random Number Generator
 */
struct RNG
{
    const std::string _name;

    RNG(std::string const &name);
    virtual ~RNG();
    /**
     * Sets the state of the RNG given a single integer value. This
     * may be used to set the state using, for example, the current time
     * stamp.
     *
     * @param seed Unsigned integer from which the initial state is generated.
     * The RNG must be able to construct a valid initial state from any
     * unsigned integer
     */
    virtual void init(unsigned int seed) = 0;
    /**
     * Returns the internal state of the RNG as a vector of integers
     *
     * @param state Vector of integers that is overwritten with the RNG state.
     */
    virtual void getState(std::vector<int> &state) const = 0;
    /**
     * Sets the internal state of the RNG given a vector of integers
     *
     * @param state Vector of integers describing the internal state of the RNG.
     * Only vectors previously derived from a call to getState should be used
     * as input to setState. Other values are not guaranteed to give a valid
     * state.
     */
    virtual bool setState(std::vector<int> const &state) = 0;
    /**
     * Generates a random value with a uniform distribution on (0,1)
     */
    virtual double uniform() = 0;
    /**
     * Generates a standard normal random value
     */
    virtual double normal() = 0;
    /**
     * Generates are andom value with an exponential distribution
     */
    virtual double exponential() = 0;
    /**
     * This static utility function may be used by an RNG object to coerce
     * values in the range [0,1] to the open range (0,1)
     */
    static double fixup(double x);
    /**
     * Returns the name of the RNG
     */
    std::string const &name() const;

    static double uniform(void *);
    static double normal(void *);
    static double exponential(void *);
};

} /* namespace jags */

#endif /* RNG_H_ */
