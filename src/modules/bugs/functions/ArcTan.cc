#include <config.h>
#include "ArcTan.h"

#include <cmath>

using std::vector;
using std::atan;
using std::string;

namespace jags {
namespace bugs {

    ArcTan::ArcTan ()
	: ScalarFunction ("arctan", 1)
    {
    }

    string ArcTan::alias() const
    {
	return "atan";
    }
    
    double ArcTan::evaluate(vector<double const *> const &args) const
    {
	return atan(*args[0]);
    }

    bool ArcTan::isDifferentiable(unsigned long i) const
    {
	return i == 0;
    }
    
    double ArcTan::gradient(vector<double const *> const &args,
			    unsigned long i) const
    {
	double z = *args[0];
	return 1/(1 + z*z);
    }

}}
