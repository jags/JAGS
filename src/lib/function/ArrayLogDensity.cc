#include <config.h>
#include <distribution/ArrayDist.h>
#include <function/ArrayLogDensity.h>
#include <util/dim.h>

using std::vector;
using std::string;

namespace jags {

    ArrayLogDensity::ArrayLogDensity(ArrayDist const *dist)
	: ArrayFunction(string("logdensity.") + dist->name().substr(1), 
			 dist->npar() + 1),
	  _dist(dist)
    {}
    
    vector<unsigned long> ArrayLogDensity::dim(
	vector<vector<unsigned long> > const &,
	vector<double const *> const &) const
    {
	return vector<unsigned long>(1, 1);
    }

    void ArrayLogDensity::evaluate(
	double *value,
	vector<double const *> const &args,
	vector<vector<unsigned long> > const &dims) const
    {
	unsigned long npar = _dist->npar();

	vector<double const *> dargs(npar);
	vector<vector<unsigned long> > ddims(npar);
	for (unsigned long i = 0; i < npar; ++i) {
	    dargs[i] = args[i+1];
	    ddims[i] = dims[i+1];
	}

	value[0] = _dist->logDensity(args[0], PDF_FULL, dargs, ddims);
    }

    

    bool
    ArrayLogDensity::checkParameterDim(vector<vector<unsigned long> > const &dim)
	const
    {
	unsigned long npar = _dist->npar();

	vector<vector<unsigned long> > ddim(npar);
	for (unsigned long i = 0; i < npar; ++i) {
	    ddim[i] = dim[i+1];
	}

	if (!_dist->checkParameterDim(ddim)) return false;
	if (dim[0] != _dist->dim(ddim)) return false;

	return true;
    }
    
    bool 
    ArrayLogDensity::checkParameterValue(vector<double const *> const &args,
					  vector<vector<unsigned long> > const &dims) 
	const
    {
	//We have to include discreteness check here as there is
	//no equivalent of checkParameterDiscrete for Functions.

	unsigned long npar = _dist->npar();

	vector<bool> mask(npar);
	for (unsigned long i = 0; i < npar; ++i) {
	    double p = *args[i + 1];
	    mask[i] = (p == static_cast<long>(p));
	}
	if (!_dist->checkParameterDiscrete(mask)) return false;

	if (_dist->isDiscreteValued(mask)) {
	    if (*args[0] != static_cast<long>(*args[0])) {
		return false;
	    }
	}

	vector<double const *> dargs(npar);
	vector<vector<unsigned long> > ddims(npar);
	for (unsigned long i = 0; i < npar; ++i) {
	    dargs[i] = args[i+1];
	    ddims[i] = dims[i+1];
	}
	return _dist->checkParameterValue(dargs, ddims);
    }

}
