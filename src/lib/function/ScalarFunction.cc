#include <config.h>
#include <function/ScalarFunction.h>
#include <util/dim.h>

#include <algorithm>

using std::vector;
using std::string;
using std::find_if;

namespace jags {

ScalarFunction::ScalarFunction (string const &name, unsigned long npar)
  : Function (name, npar)
{
}

bool 
ScalarFunction::checkParameterValue(vector<double const *> const &) const
{
    return true;
}

bool ScalarFunction::isPower(vector<bool> const &mask,
			     vector<bool> const &) const
{
    unsigned long nmask = 0;
    for (unsigned long i = 0; i < mask.size(); ++i) {
	nmask += mask[i];
    }
    
    if (nmask > 1)
	return false;
    else
	return isScale(mask, vector<bool>());
    
}

} //namespace jags
