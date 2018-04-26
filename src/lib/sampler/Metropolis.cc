#include <config.h>
#include <sampler/Metropolis.h>
#include <rng/RNG.h>

#include <algorithm>
#include <stdexcept>

using std::logic_error;
using std::vector;
using std::copy;
using std::min;

namespace jags {

Metropolis::Metropolis(vector<double> const &value)
    : _last_value(value), _adapt(true)
{
}

Metropolis::~Metropolis()
{
}

bool Metropolis::accept(RNG *rng, double prob)
{
    bool accept = rng->uniform() <= prob;
    if (accept) {
	//Store current value as last accepted value
	getValue(_last_value);
    }
    else {
	//Revert to last accepted value
	setValue(_last_value);
    }
    if (_adapt) {
	rescale(min(prob, 1.0));
    }
    return accept;
}

void Metropolis::adaptOff()
{
    _adapt = false;
}

bool Metropolis::isAdaptive() const
{
    return true;
}

unsigned long Metropolis::length() const
{
    return _last_value.size();
}

} //namespace jags
