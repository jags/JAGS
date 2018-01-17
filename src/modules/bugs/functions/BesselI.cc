#include <config.h>

#include "BesselI.h"

#include <JRmath.h>

using std::vector;

namespace jags {
    namespace bugs {

	BesselI::BesselI ()
	    : ScalarFunction ("besselI", 3)
	{
	}

	double BesselI::evaluate(vector<double const *> const &args) const
	{
	    unsigned int expo = 1 + (*args[2] != 0);
	    return bessel_i(*args[0], *args[1], expo);
	}

	bool BesselI::checkParameterValue(vector <double const *> const &args)
	{
	    return *args[0] >= 0;
	}

    }
}
