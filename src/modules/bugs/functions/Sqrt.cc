#include <config.h>
#include "Sqrt.h"

#include <cmath>

using std::vector;
using std::sqrt;

namespace jags {
namespace bugs {

    Sqrt::Sqrt ():ScalarFunction ("sqrt", 1)
    {
    }

    double Sqrt::evaluate(vector<double const *> const &args) const
    {
	return sqrt(*args[0]);
    }

    bool Sqrt::checkParameterValue(vector<double const *> const &args) const
    {
	return *args[0] >= 0;
    }

    bool Sqrt::isPower(vector<bool> const &, vector<bool> const &) const
    {
        return true;
    }

    bool Sqrt::isDifferentiable(unsigned long i) const
    {
	return i == 0;
    }
    
    double Sqrt::gradient(vector<double const *> const &args,
			 unsigned long i) const
    {
	return 1/(2*sqrt(*args[0]));
    }

}}
