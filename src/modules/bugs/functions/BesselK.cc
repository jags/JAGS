#include <config.h>

#include "BesselK.h"

#include <JRmath.h>

using std::vector;

namespace jags {
    namespace bugs {

	BesselK::BesselK ()
	    : ScalarFunction ("BesselK", 3)
	{
	}

	double BesselK::evaluate(vector<double const *> const &args) const
	{
	    unsigned int expo = 1 + (*args[2] != 0);
	    return bessel_k(*args[0], *args[1], expo);
	}

	bool BesselK::checkParameterValue(vector <double const *> const &args)
	{
	    return *args[0] >= 0;
	}

    }
}
