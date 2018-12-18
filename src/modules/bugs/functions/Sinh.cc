#include <config.h>
#include "Sinh.h"

#include <cmath>

using std::vector;
using std::sinh;

namespace jags {
namespace bugs {

    Sinh::Sinh ()
	: ScalarFunction ("sinh", 1)
    {
    }

    double Sinh::evaluate(vector<double const *> const &args) const
    {
	return sinh(*args[0]);
    }

    bool Sinh::isDifferentiable(unsigned long i) const
    {
	return i == 0;
    }
    
    double Sinh::gradient(vector<double const *> const &args,
			  unsigned long i) const
    {
	return cosh(*args[0]);
    }

}}
