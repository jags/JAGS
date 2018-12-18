#include <config.h>
#include "Tan.h"

#include <cmath>

using std::vector;
using std::tan;

namespace jags {
namespace bugs {


    Tan::Tan ()
	: ScalarFunction ("tan", 1)
    {
    }

    double Tan::evaluate(vector<double const *> const &args) const
    {
	return tan(*args[0]);
    }

    bool Tan::isDifferentiable(unsigned long i) const
    {
	return i == 0;
    }
    
    double Tan::gradient(vector<double const *> const &args,
			 unsigned long i) const
    {
	double y = tan(*args[0]);
	return 1 + y*y;
    }

}}
