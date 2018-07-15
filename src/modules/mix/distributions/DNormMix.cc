#include <config.h>
#include "DNormMix.h"
#include <util/dim.h>
#include <util/nainf.h>

#include <cmath>

#include <JRmath.h>

using std::vector;

#define MU(par) (par[0])
#define TAU(par) (par[1])
#define PROB(par) (par[2])
#define NCAT(lengths) (lengths[0])

namespace jags {
namespace mix {

    DNormMix::DNormMix()
	: VectorDist("dnormmix", 3)
    {}

    bool 
    DNormMix::checkParameterLength(vector<unsigned long> const &lengths) const
    {
	if (lengths[0] == 1)
	    return false; // Must be a mixture
   
	// Parameter lengths must match 
	return (lengths[0] == lengths[1]) && (lengths[0] == lengths[2]);
    }

    bool 
    DNormMix::checkParameterValue(vector<double const *> const &par,
				  vector<unsigned long> const &lengths) const
    {
	for (unsigned int i = 0; i < NCAT(lengths); ++i) {
	    if (TAU(par)[i] <= 0)
		return false;
	    if (PROB(par)[i] <= 0) 
		return false;
	}
	return true;
    }

    double DNormMix::logDensity(double const *x,
				PDFType type,
				vector<double const *> const &par,
				vector<unsigned long> const &lengths) const
    {
	double density = 0.0;
	double psum = 0.0;
	for (unsigned int i = 0; i < NCAT(lengths); ++i) {
	    density += PROB(par)[i] * dnorm(*x, MU(par)[i], 
					    1/sqrt(TAU(par)[i]), 0);
	    psum += PROB(par)[i];
	}
	return log(density) - log(psum);
    }

    void 
    DNormMix::randomSample(double *x,
			   vector<double const *> const &par, 
			   vector<unsigned long> const &lengths, RNG *rng) const
    {
	unsigned long Ncat = NCAT(lengths);

	// Rescale probability parameter
	double sump = 0;
	for (unsigned int i = 0; i < Ncat; ++i) {
	  sump += PROB(par)[i];
	}
	double p_rand = runif(0, 1, rng) * sump;
	sump = 0;

	// Select mixture component (r)
	unsigned long r = Ncat - 1;
	for (unsigned long i = 0; i < Ncat - 1; ++i) {
	    sump += PROB(par)[i];
	    if (sump > p_rand) {
		r = i;
		break;
	    }
	}

	// Now sample from conditional distribution of component r
	double ans = rnorm(MU(par)[r], 1/sqrt(TAU(par)[r]), rng);
	*x = ans;
    }

    bool DNormMix::isSupportFixed(vector<bool> const &fixmask) const
    {
	return true;
    }

    void 
    DNormMix::support(double *lower, double *upper,
		      vector<double const *> const &parameters,
		      vector<unsigned long> const &lengths) const
    {
	*lower = JAGS_NEGINF;
	*upper = JAGS_POSINF;
    }

    unsigned long DNormMix::length(vector<unsigned long> const &parlengths) const
    {
	return 1;
    }

}}
