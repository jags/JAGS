#include "BaseRNGFactory.h"
#include "WichmannHillRNG.h"
#include "MarsagliaRNG.h"
#include "SuperDuperRNG.h"
#include "MersenneTwisterRNG.h"

#include <ctime>

using std::vector;

#define DEFAULT_NORM_KIND KINDERMAN_RAMAGE

BaseRNGFactory::BaseRNGFactory()
    : _index(0)
{
}

BaseRNGFactory::~BaseRNGFactory()
{
  for (unsigned int i = 0; i < _rngvec.size(); ++i) {
    delete _rngvec[i];
  }
}

RNG * BaseRNGFactory::makeRNG()
{
    unsigned int seed = static_cast<unsigned int>(time(NULL));

    RNG *rng = 0;

    switch(_index++) {
    case 0:
	rng =  new WichmannHillRNG(seed, DEFAULT_NORM_KIND);
	break;
    case 1:
	rng = new MarsagliaRNG(seed, DEFAULT_NORM_KIND);
	break;
    case 2:
	rng = new SuperDuperRNG(seed, DEFAULT_NORM_KIND);
	break;
    case 3:
        rng = new MersenneTwisterRNG(seed, DEFAULT_NORM_KIND);
	break;
    default:
	return 0;
    }

    // Store generated RNG for memory management
    _rngvec.push_back(rng);
    return rng;
}

RNG * BaseRNGFactory::makeRNG(std::string const &name)
{
    unsigned int seed = static_cast<unsigned int>(time(NULL));

    RNG *rng = 0;
    if (name == "base::Wichmann-Hill")
	rng = new WichmannHillRNG(seed, DEFAULT_NORM_KIND);
    else if (name == "base::Marsaglia-Multicarry")
	rng = new MarsagliaRNG(seed, DEFAULT_NORM_KIND);
    else if (name == "base::Super-Duper")
	rng = new SuperDuperRNG(seed, DEFAULT_NORM_KIND);
    else if (name == "base::Mersenne-Twister")
        rng = new MersenneTwisterRNG(seed, DEFAULT_NORM_KIND);
    else
	return 0;

    _rngvec.push_back(rng);
    return rng;
}
