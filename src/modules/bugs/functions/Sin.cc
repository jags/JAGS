#include <config.h>
#include "Sin.h"

#include <cmath>

using std::vector;
using std::sin;

namespace jags {
namespace bugs {


    Sin::Sin ()
	: ScalarFunction ("sin", 1)
    {
    }

    double Sin::evaluate(vector<double const *> const &args) const
    {
	return sin(*args[0]);
    }

    bool Sin::isDifferentiable(unsigned long i) const
    {
	return i == 0;
    }
    
    double Sin::gradient(vector<double const *> const &args,
			 unsigned long i) const
    {
	return cos(*args[0]);
    }

}}
