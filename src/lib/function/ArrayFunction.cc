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
ArrayFunction::checkParameterValue(vector<double const *> const &args,
				   vector<vector<unsigned long> > const &dims) 
    const
{
    return true;
}

} //namespace jags
