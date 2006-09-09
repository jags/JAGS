#ifndef RNG_FACTORY_H_
#define RNG_FACTORY_H_

#include <string>

class RNG;

/**
 * @short Factory object for Random Number Generators
 */
class RNGFactory
{
 public:
  /**
   * Destructor. An RNGFactory retains ownership of the RNG objects
   * it generates, and should delete them when the destructor is called.
   */
  virtual ~RNGFactory() {};
  /**
   * Returns a newly allocated RNG object.
   *
   * Successive calls should return RNG objects that generate
   * completely independent random number streams.  An RNGFactory may
   * have the capacity to generate only a limited number of
   * independent RNG objects.  When this capacity is exhausted, the
   * RNGFactory will return a NULL pointer.
   */
  virtual RNG * makeRNG() = 0;
  /**
   * Returns a newly allocated RNG object.  
   *
   * This function can be repeatedly called with the same name
   * argument. There is no guarantee that RNG objects created in this
   * way will generate independent streams.
   */
  virtual RNG * makeRNG(std::string const &name) = 0;
};

#endif /* RNG_FACTORY_H_ */
