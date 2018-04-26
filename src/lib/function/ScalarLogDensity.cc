#include <config.h>
#include <distribution/ScalarDist.h>
#include <function/ScalarLogDensity.h>

using std::vector;
using std::string;

namespace jags {

    ScalarLogDensity::ScalarLogDensity(ScalarDist const *dist)
	: ScalarFunction(string("logdensity.") + dist->name().substr(1), 
			 dist->npar() + 1),
	  _dist(dist)
    {}
    
    double ScalarLogDensity::evaluate(vector<double const *> const &args) const
    {
	vector<double const *> param(_dist->npar());
	for (unsigned long i = 0; i < param.size(); ++i) {
	    param[i] = args[i+1];
	}

	return _dist->logDensity(*args[0], PDF_FULL, param, 0, 0);
    }

    bool 
    ScalarLogDensity::checkParameterValue(vector<double const *> const &args) 
	const
    {
	//We have to include discreteness check here as there is
	//no equivalent of checkParameterDiscrete for Functions.

	vector<bool> mask(_dist->npar());
	for (unsigned long i = 0; i < mask.size(); ++i) {
	    double p = *args[i + 1];
	    mask[i] = (p == static_cast<int>(p));
	}
	if (!_dist->checkParameterDiscrete(mask)) return false;

	if (_dist->isDiscreteValued(mask)) {
	    if (*args[0] != static_cast<int>(*args[0])) {
		return false;
	    }
	}

	vector<double const *> param(_dist->npar());
	for (unsigned long i = 0; i < param.size(); ++i) {
	    param[i] = args[i+1];
	}
	return _dist->checkParameterValue(param);
    }

}
