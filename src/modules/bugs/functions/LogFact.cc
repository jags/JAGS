#include <config.h>
#include "LogFact.h"

#include <JRmath.h>

using std::vector;

namespace jags {
namespace bugs {

    LogFact::LogFact ()
	: ScalarFunction ("logfact", 1)
    {
    }

    double LogFact::evaluate(vector<double const *> const &args) const
    {
	return lgammafn(*args[0] + 1);
    }

    bool LogFact::checkParameterValue(vector<double const *> const &args) const
    {
	return *args[0] > -1;
    }

    bool LogFact::isDifferentiable(unsigned long i) const
    {
	return i == 0;
    }
    
    double LogFact::gradient(vector<double const *> const &args,
			     unsigned long i) const
    {
	return digamma(*args[0] + 1);
    }

}}
