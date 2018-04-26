#include <config.h>
#include "DIntervalFunc.h"
#include <util/dim.h>

using std::vector;

#define T(args) (*args[0])
#define CUTPOINTS(args) (args[1])

namespace jags {
namespace bugs {

    DIntervalFunc::DIntervalFunc () : ScalarVectorFunction ("dinterval", 2)
    {
    }

    double DIntervalFunc::scalarEval(vector<double const *> const &args,
				     vector<unsigned long> const &lengths) const
    {
	unsigned long ncut = lengths[1];
	double t = T(args);
	for (unsigned long i = 0; i < ncut; ++i) {
	    if (t <= CUTPOINTS(args)[i])
		return i;
	}
	return ncut;
    }
    
    bool DIntervalFunc::checkParameterLength (vector<unsigned long> const &args)
	const
    {
	return args[0] == 1 && args[1] >= 1;
    }
    
    bool DIntervalFunc::isDiscreteValued(vector<bool> const &mask) const
    {
	return true;
    }
    
}}
