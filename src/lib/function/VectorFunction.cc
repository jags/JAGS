#include <config.h>
#include <function/VectorFunction.h>

using std::vector;
using std::string;

namespace jags {

VectorFunction::VectorFunction (string const &name, unsigned long npar)
    : Function(name, npar)
{
}

bool 
VectorFunction::checkParameterValue(vector<double const *> const &,
				    vector<unsigned long> const &) const
{
    return true;
}

bool 
VectorFunction::checkParameterLength(vector<unsigned long> const &arglen) const
{
    for (unsigned long i = 0; i < arglen.size(); ++i) {
	if (arglen[i] == 0) return false;
    }
    return true;
}

    void
    VectorFunction::gradient(double *grad, vector<double const *> const &args,
			     vector<unsigned long> const &dims,
			     unsigned long i) const
    {
    }

} //namespace jags
