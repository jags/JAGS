#include <config.h>
#include "Round.h"

#include <cmath>

using std::vector;
using std::floor;

namespace jags {
namespace bugs {

    Round::Round ():ScalarFunction ("round", 1)
    {
    }

    double Round::evaluate(vector<double const *> const &args) const
    {
	/* The C99 round function does not exist in the current C++ standard */
	//FIXME: C++11
	return floor(*args[0] + 0.5);
    }

    bool Round::isDiscreteValued(vector<bool> const &) const
    {
	return true;
    }

}}
