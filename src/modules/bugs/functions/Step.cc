#include <config.h>
#include "Step.h"

using std::vector;

namespace jags {
namespace bugs {

    Step::Step ()
	: ScalarFunction ("step", 1)
    {
    }

    double
    Step::evaluate(vector<double const *> const &args) const
    {
	return *args[0] >= 0 ? 1 : 0;
    }

    bool Step::isDiscreteValued(vector<bool> const &) const
    {
	return true;
    }

}}
