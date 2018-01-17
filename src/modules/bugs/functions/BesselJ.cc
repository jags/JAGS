#include <config.h>

#include "BesselJ.h"

#include <JRmath.h>

using std::vector;

namespace jags {
    namespace bugs {

	BesselJ::BesselJ ()
	    : ScalarFunction ("BesselJ", 2)
	{
	}

	double BesselJ::evaluate(vector<double const *> const &args) const
	{
	    return bessel_j(*args[0], *args[1]);
	}

	bool BesselJ::checkParameterValue(vector <double const *> const &args)
	{
	    return *args[0] >= 0;
	}

    }
}
