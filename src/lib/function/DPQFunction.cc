#include <config.h>
#include <function/DPQFunction.h>
#include <distribution/RScalarDist.h>

using std::vector;
using std::string;

namespace jags {

    DPQFunction::DPQFunction(string const &name, RScalarDist const *dist) 
	: ScalarFunction(name, dist->npar() + 1), _dist(dist)
    {}
    
    RScalarDist const *DPQFunction::dist() const
    {
	return _dist;
    }
    
    bool DPQFunction::checkArgs(vector<double const *> const &args) const
    {
	vector<double const *> param(_dist->npar());
	for (unsigned long i = 0; i < param.size(); ++i) {
	    param[i] = args[i+1];
	}
	
	return _dist->checkParameterValue(param);
    }

}
