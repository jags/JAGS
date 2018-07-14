#include <config.h>
#include "Equals.h"

#include <cmath>
#include <cfloat>

using std::vector;
using std::fabs;

namespace jags {
namespace bugs {

    Equals::Equals () : ScalarFunction ("equals",2)
    {
    }

    double Equals::evaluate(vector<double const *> const &args) const
    {
	return fabs(*args[0] - *args[1]) < 16 *DBL_EPSILON;
    }

    bool Equals::isDiscreteValued(vector<bool> const &) const
    {
	return true;
    }

}}
