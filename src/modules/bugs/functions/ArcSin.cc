#include <config.h>
#include "ArcSin.h"
#include <util/nainf.h>

#include <cmath>

using std::vector;
using std::asin;
using std::string;

namespace jags {
namespace bugs {

    ArcSin::ArcSin ()
	: ScalarFunction ("arcsin", 1)
    {
    }

    string ArcSin::alias() const
    {
	return "asin";
    }
    
    double ArcSin::evaluate(vector<double const *> const &args) const
    {
	return asin(*args[0]);
    }

    bool ArcSin::checkParameterValue(vector<double const *> const &args) const
    {
	return *args[0] >= -1 && *args[0] <= 1;
    }

    bool ArcSin::isDifferentiable(unsigned long i) const
    {
	return i == 0;
    }
    
    double ArcSin::gradient(vector<double const *> const &args,
			    unsigned long i) const
    {
	double z = *args[0];
	return 1/sqrt(1 - z*z);
    }

}}
