#include <config.h>
#include <distribution/DistScalarRmath.h>
#include <rng/RNG.h>
#include <util/nainf.h>
#include <util/dim.h>

#include <stdexcept>
#include <cmath>
#include <algorithm>

using std::string;
using std::vector;
using std::length_error;
using std::logic_error;
using std::log;
using std::min;
using std::max;

double DistScalarRmath::calPlower(double lower, 
				  vector<double const*> const &parameters) const
{
    if (isDiscreteValued()) {
	return p(lower - 1, parameters, true, false);
    }
    else {
	return p(lower, parameters, true, false);
    }
}

double DistScalarRmath::calPupper(double upper,
				  vector<double const*> const &parameters) const
{
    return p(upper, parameters, true, false);
}


DistScalarRmath::DistScalarRmath(string const &name, unsigned int npar, 
				 Support support, bool canbound, bool discrete)
  
    : DistScalar(name, npar, support, canbound, discrete),  _support(support)
{
}

double 
DistScalarRmath::typicalScalar(vector<double const *> const &parameters,
			       double const *lower, double const *upper) const
{
    double llimit = JAGS_NEGINF, ulimit = JAGS_POSINF;
    support(&llimit, &ulimit, parameters);
    double plower = 0, pupper = 1;
    
    if (lower) {
	llimit = max(llimit, *lower);
	plower = calPlower(llimit, parameters);
    }

    if (upper) {
	ulimit = min(ulimit, *upper);
	pupper = calPupper(ulimit, parameters);
    }
    
    double pmed = (plower + pupper)/2;
    double med = q(pmed, parameters, true, false);	

    //Calculate the log densities
    double dllimit = d(llimit, parameters, true);
    double dulimit = d(ulimit, parameters, true);
    double dmed = d(med, parameters, true);

    //Pick the median if it has the highest density, otherwise pick
    //a point near to (but not on) the boundary
    if (dmed >= dllimit && dmed >= dulimit) {
	return med;
    }
    else if (dulimit > dllimit) {
	return q(0.1 * plower + 0.9 * pupper, parameters, true, false);
    }
    else {
	return q(0.9 * plower + 0.1 * pupper, parameters, true, false);
    }
}

double 
DistScalarRmath::scalarLogLikelihood(double x,
				     vector<double const *> const &parameters,
				     double const *lower, double const *upper) 
  const
{
    double loglik =  d(x, parameters, true);

    if (lower || upper) {

	double plower = 0, pupper = 1;

	if (lower) {

	    if (x < *lower)
		return JAGS_NEGINF;

	    plower = calPlower(*lower, parameters);
	}

	if (upper) {
	    
	    if (x > *upper)
		return JAGS_NEGINF;

	    pupper = calPupper(*upper, parameters);
	}

	loglik -= log(pupper - plower);
    }

    return loglik;
}


double 
DistScalarRmath::scalarRandomSample(vector<double const *> const &parameters,
				    double const *lower, double const *upper,
				    RNG *rng) const
{
    if (lower || upper) {

	double plower = 0, pupper = 1;
	if (lower) {
	    plower = calPlower(*lower, parameters);
	}
	if (upper) {
	    pupper = calPupper(*upper, parameters);
	}
	
	double u = plower + rng->uniform() * (pupper - plower);
	return q(u, parameters, true, false);
    }
    else {
	return r(parameters, rng);
    }

}
