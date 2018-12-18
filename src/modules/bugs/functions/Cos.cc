#include <config.h>
#include "Cos.h"

#include <cmath>

using std::vector;
using std::cos;

namespace jags {
namespace bugs {

    Cos::Cos ()
	: ScalarFunction ("cos", 1)
    {
    }

    double Cos::evaluate(vector<double const *> const &args) const
    {
	return cos(*args[0]);
    }

    bool Cos::isDifferentiable(unsigned long i) const
    {
	return i == 0;
    }
    
    double Cos::gradient(vector<double const *> const &args,
			 unsigned long index) const
    {
	return -sin(*args[0]);
    }

}}
