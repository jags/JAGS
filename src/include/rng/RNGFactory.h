#ifndef RNG_FACTORY_H_
#define RNG_FACTORY_H_

#include <string>
#include <vector>

namespace jags {

struct RNG;

/**
 * @short Factory for RNG objects
 */
class RNGFactory
{
    
  public:
    /**
     * Destructor. An RNGFactory retains ownership of the RNG objects
     * it generates, and should delete them when the destructor is called.
     */
    virtual ~RNGFactory();
    /**
     * Sets the random seed of the RNG factory so that a reproducible
     * sequence of RNGs can be produced.
     *
     * @param seed Seed that uniquely determines the sequence of RNGs 
     * produced by subsequent calls to makeRNGs.
     */
    virtual void setSeed(unsigned int seed) = 0;
    /**
     * Returns a vector of newly allocated RNG objects.
     *
     * @param n Number of RNGs requested. Note that an RNG factory may have the
     * capacity to make only m < n independent samplers (where m may be zero).
     * In this case it should return a vector of length m.
     */
    virtual std::vector<RNG *> makeRNGs(unsigned int n) = 0;
    /**
     * Returns a newly allocated RNG object.  
     *
     * This function can be repeatedly called with the same name
     * argument. There is no guarantee that RNG objects created in this
     * way will generate independent streams.
     */
    virtual RNG * makeRNG(std::string const &name) = 0;
    /**
     * Returns the name of the RNG factory
     */
    virtual std::string name() const = 0;
};

} /* namespace jags */

#endif /* RNG_FACTORY_H_ */
