#include <config.h>
#include "ArcCos.h"

#include <cmath>

using std::vector;
using std::acos;
using std::string;

namespace jags {
namespace bugs {

    ArcCos::ArcCos ()
	: ScalarFunction ("arccos", 1)
    {
    }

    string ArcCos::alias() const
    {
	return "acos";
    }
    
    double ArcCos::evaluate(vector<double const *> const &args) const
    {
	return acos(*args[0]);
    }

    bool ArcCos::checkParameterValue(vector<double const *> const &args) const
    {
	return *args[0] >= -1 && *args[0] <= 1;
    }

    bool ArcCos::isDifferentiable(unsigned long i)
    {
	return i == 0;
    }
    
    double ArcCos::gradient(vector<double const *> const &args,
			    unsigned long i) const
    {
	double z = *args[0];
	return -1/sqrt(1 - z*z);
    }

}}
