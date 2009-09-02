
#include <config.h>
#include "Trunc.h"

#include <cmath>

using std::vector;
using std::floor;

namespace bugs {

    Trunc::Trunc ():ScalarFunc ("trunc", 1)
    {
    }

    double Trunc::evaluateScalar(vector<double const *> const &args) const
    {
	/* The C99 trunc function does not exist in the current C++ standard. */
	return args[0][0] >= 0 ? floor(args[0][0]) : -floor(-args[0][0]);
    }

    bool Trunc::isDiscreteValued(vector<bool> const &mask) const
    {
	return true;
    }
}
