#include <config.h>

#include "BesselY.h"

#include <JRmath.h>

using std::vector;

namespace jags {
    namespace bugs {

	BesselY::BesselY ()
	    : ScalarFunction ("BesselY", 2)
	{
	}

	double BesselY::evaluate(vector<double const *> const &args) const
	{
	    return bessel_y(*args[0], *args[1]);
	}

	bool BesselY::checkParameterValue(vector <double const *> const &args)
	{
	    return *args[0] >= 0;
	}

    }
}
