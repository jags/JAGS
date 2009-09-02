#include <config.h>
#include "Round.h"

#include <cmath>

using std::vector;
using std::floor;

namespace bugs {

    Round::Round ():ScalarFunc ("round", 1)
    {
    }

    double Round::evaluateScalar(vector<double const *> const &args) const
    {
	/* The C99 round function does not exist in the current C++ standard */
	return floor(*args[0] + 0.5);
    }

    bool Round::isDiscreteValued(vector<bool> const &mask) const
    {
	return true;
    }
}
