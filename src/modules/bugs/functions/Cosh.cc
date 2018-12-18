#include <config.h>
#include "Cosh.h"

#include <cmath>

using std::vector;
using std::cosh;

namespace jags {
namespace bugs {

    Cosh::Cosh ()
	: ScalarFunction ("cosh", 1)
    {
    }

    double Cosh::evaluate(vector<double const *> const &args) const
    {
	return cosh(*args[0]);
    }

    bool Cosh::isDifferentiable(unsigned long i) const
    {
	return i == 0;
    }
    
    double Cosh::gradient(vector<double const *> const &args,
			  unsigned long i) const
    {
	return sinh(*args[0]);
    }

}}
