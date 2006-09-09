#ifndef BASE_RNG_FACTORY_H_
#define BASE_RNG_FACTORY_H_

#include <rng/RNGFactory.h>

#include <vector>

class RNG;

/**
 * @short Factory object for Random Number Generators
 */
class BaseRNGFactory : public RNGFactory
{
    unsigned int _index;
    std::vector<RNG*> _rngvec;
public:
    BaseRNGFactory();
    ~BaseRNGFactory();
    RNG * makeRNG();
    RNG * makeRNG(std::string const &name);
};

#endif /* BASE_RNG_FACTORY_H_ */
