#include <config.h>
#include "Tanh.h"

#include <cmath>

using std::vector;
using std::tanh;

namespace jags {
namespace bugs {

    Tanh::Tanh ()
	: ScalarFunction ("tanh", 1)
    {
    }

    double Tanh::evaluate(vector<double const *> const &args) const
    {
	return tanh(*args[0]);
    }

    bool Tanh::isDifferentiable(unsigned long i) const
    {
	return i == 0;
    }
    
    bool Tanh::gradient(double &grad, vector<double const *> const &args,
			unsigned long i) const
    {
	double y = cosh(*args[0]);
	grad = 1/(y*y);
	return true;
    }

}}
