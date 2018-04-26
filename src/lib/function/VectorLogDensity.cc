#include <config.h>
#include <distribution/VectorDist.h>
#include <function/VectorLogDensity.h>

using std::vector;
using std::string;

namespace jags {

    VectorLogDensity::VectorLogDensity(VectorDist const *dist)
	: VectorFunction(string("logdensity.") + dist->name().substr(1), 
			 dist->npar() + 1),
	  _dist(dist)
    {}

    unsigned long 
    VectorLogDensity::length(vector<unsigned long> const &lengths,
			     vector<double const *> const &values) const
    {
	return 1;
    }
    
    void 
    VectorLogDensity::evaluate(double *value,
			       vector<double const *> const &args,
			       vector<unsigned long> const &lengths) const
    {
	unsigned long npar = _dist->npar();

	vector<double const *> dargs(npar);
	vector<unsigned long > dlengths(npar);
	for (unsigned long i = 0; i < npar; ++i) {
	    dargs[i] = args[i+1];
	    dlengths[i] = lengths[i+1];
	}

	value[0] = _dist->logDensity(args[0], lengths[0], PDF_FULL, 
				     dargs, dlengths, 0, 0);
    }

    bool
    VectorLogDensity::checkParameterLength(vector<unsigned long> const &lengths)
	const
    {
	unsigned long npar = _dist->npar();

	vector<unsigned long> dlengths(npar);
	for (unsigned long i = 0; i < npar; ++i) {
	    dlengths[i] = lengths[i+1];
	}

	if (!_dist->checkParameterLength(dlengths)) return false;
	if (lengths[0] != _dist->length(dlengths)) return false;

	return true;
    }

    bool 
    VectorLogDensity::checkParameterValue(vector<double const *> const &args,
					  vector<unsigned long> const &lengths) 
	const
    {
	//We have to include discreteness check here as there is
	//no equivalent of checkParameterDiscrete for Functions.

	unsigned long npar = _dist->npar();

	vector<bool> mask(npar);
	for (unsigned long i = 0; i < npar; ++i) {
	    double p = *args[i + 1];
	    mask[i] = (p == static_cast<int>(p));
	}
	if (!_dist->checkParameterDiscrete(mask)) return false;

	if (_dist->isDiscreteValued(mask)) {
	    if (*args[0] != static_cast<int>(*args[0])) {
		return false;
	    }
	}

	vector<double const *> dargs(npar);
	vector<unsigned long> dlengths(npar);
	for (unsigned long i = 0; i < npar; ++i) {
	    dargs[i] = args[i+1];
	    dlengths[i] = lengths[i+1];
	}
	return _dist->checkParameterValue(dargs, dlengths);
    }

}
