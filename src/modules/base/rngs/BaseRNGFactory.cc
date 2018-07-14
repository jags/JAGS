#include "BaseRNGFactory.h"
#include "WichmannHillRNG.h"
#include "MarsagliaRNG.h"
#include "SuperDuperRNG.h"
#include "MersenneTwisterRNG.h"

#include <ctime>
#include <climits>

using std::vector;
using std::time;
using std::string;

#define DEFAULT_NORM_KIND KINDERMAN_RAMAGE

namespace jags {
namespace base {

    BaseRNGFactory::BaseRNGFactory()
	: _index(0), _seed(static_cast<unsigned int>(time(nullptr)))
    {
    }

    BaseRNGFactory::~BaseRNGFactory()
    {
	for (unsigned int i = 0; i < _rngvec.size(); ++i) {
	    delete _rngvec[i];
	}
    }

    void BaseRNGFactory::setSeed(unsigned int seed)
    {
	_seed = seed;
    }

    vector<RNG *> BaseRNGFactory::makeRNGs(unsigned int n)
    {
	vector<RNG *> ans;
	for (unsigned int i = 0; i < n; i++) {
	    RNG *rng = nullptr;
	    switch(_index) {
	    case 0:
		rng =  new WichmannHillRNG(_seed, DEFAULT_NORM_KIND);
		break;
	    case 1:
		rng = new MarsagliaRNG(_seed, DEFAULT_NORM_KIND);
		break;
	    case 2:
		rng = new SuperDuperRNG(_seed, DEFAULT_NORM_KIND);
		break;
	    case 3:
		rng = new MersenneTwisterRNG(_seed, DEFAULT_NORM_KIND);
		break;
	    default:
		break;
	    }

	    //Move onto the next generator
	    _index = (_index + 1) % 4;

	    //Reset the seed
	    _seed = static_cast<unsigned int>(rng->uniform() * UINT_MAX);
	    
	    // Store generated RNG for memory management
	    _rngvec.push_back(rng);
	    ans.push_back(rng);

	}
	return ans;
    }

    RNG * BaseRNGFactory::makeRNG(string const &name)
    {
	unsigned int seed = static_cast<unsigned int>(time(nullptr));

	RNG *rng = nullptr;
	if (name == "base::Wichmann-Hill")
	    rng = new WichmannHillRNG(seed, DEFAULT_NORM_KIND);
	else if (name == "base::Marsaglia-Multicarry")
	    rng = new MarsagliaRNG(seed, DEFAULT_NORM_KIND);
	else if (name == "base::Super-Duper")
	    rng = new SuperDuperRNG(seed, DEFAULT_NORM_KIND);
	else if (name == "base::Mersenne-Twister")
	    rng = new MersenneTwisterRNG(seed, DEFAULT_NORM_KIND);
	else
	    return nullptr;

	_rngvec.push_back(rng);
	return rng;
    }

    string BaseRNGFactory::name() const
    {
	return "base::BaseRNG";
    }

}}
