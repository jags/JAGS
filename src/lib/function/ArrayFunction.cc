#include <config.h>
#include <function/ArrayFunction.h>

using std::vector;
using std::string;

namespace jags {

ArrayFunction::ArrayFunction (string const &name, unsigned long npar)
    : Function(name, npar)
{
}

bool 
ArrayFunction::checkParameterValue(vector<double const *> const &,
				   vector<vector<unsigned long> > const &) 
    const
{
    return true;
}

} //namespace jags
