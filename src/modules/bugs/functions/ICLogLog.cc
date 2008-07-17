#include <config.h>
#include "ICLogLog.h"

#include <cmath>

using std::vector;
using std::log;
using std::exp;

namespace bugs {

    ICLogLog::ICLogLog(): InverseLinkFunc("icloglog", "cloglog")
    {
    }

    double
    ICLogLog::evaluateScalar(vector <double const *> const &args) const
    {
	return 1 - exp(-exp(*args[0]));
    }

    double ICLogLog::link(double mu) const
    {
	return log (-log (1 - mu));
    }

    double ICLogLog::grad(double eta) const
    {
	return exp(eta) * exp(-exp(eta));
    }

}
