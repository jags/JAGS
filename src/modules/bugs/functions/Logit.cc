#include <config.h>
#include "Logit.h"

#include <cmath>

using std::vector;
using std::log;

namespace jags {
namespace bugs {

    Logit::Logit ():ScalarFunction ("logit", 1)
    {
    }

    double Logit::evaluate(vector <double const *> const &args) const
    {
	double arg = *args[0];
	return log(arg) - log(1 - arg);
    }

    bool Logit::checkParameterValue (vector <double const *> const &args) const
    {
	double arg = *args[0];
	return (arg >= 0 && arg <= 1);
    }

    bool Logit::isDifferentiable(unsigned long i) const
    {
	return i == 0;
    }
    
    double Logit::gradient(vector <double const *> const &args,
			   unsigned long i) const
    {
	double x = *args[0];
	return 1/(x*(1 - x));
    }

}}
