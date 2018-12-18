#include <config.h>
#include "Abs.h"
#include <util/nainf.h>

#include <cmath>

using std::vector;
using std::fabs;

namespace jags {
namespace bugs {

    Abs::Abs ():ScalarFunction ("abs", 1)
    {
    }

    double Abs::evaluate(vector<double const *> const &args) const
    {
	return fabs(*args[0]);
    }


    bool Abs::isDiscreteValued(vector<bool> const &mask) const
    {
	return mask[0];
    }

}}
