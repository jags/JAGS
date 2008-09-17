#include <config.h>
#include "DCat.h"
#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>

#include <cfloat>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <JRmath.h>

using std::min;
using std::max;
using std::vector;
using std::max_element;
using std::logic_error;

#define PROB(par) (par[0])
#define NCAT(dims) (dims[0][0])

DCat::DCat()
    : Distribution("dcat", 1, false, true) 
{}

bool DCat::checkParameterDim(vector<vector<unsigned int> > const &dims) const
{
    return isVector(dims[0]) || isScalar(dims[0]);
}

bool DCat::checkParameterValue(vector<double const *> const &par,
			       vector<vector<unsigned int> > const &dims) const
{
    bool nz = false;
    for (unsigned int i = 0; i < NCAT(dims); i++) {
	if (PROB(par)[i] < 0.0) {
	    return false; 
	}
	if (PROB(par)[i] > 0.0) {
	    nz = true; //at least one non-zero probability
	}
    }
    return nz;
}

double DCat::logLikelihood(double const *x, unsigned int length,
			   vector<double const *> const &par,
			   vector<vector<unsigned int> > const &dims,
			   double const *lower, double const *upper) const
{
    unsigned int y = static_cast<unsigned int>(*x);
    if (y < 1 || y > NCAT(dims)) {
	return JAGS_NEGINF;
    }
    else {
	double sump = 0.0;
	for (unsigned int i = 0; i < NCAT(dims); i++) {
	    sump += PROB(par)[i];
	}
	return log(PROB(par)[y-1]) - log(sump);
    }
}

void DCat::randomSample(double *x, unsigned int length,
			vector<double const *> const &par,
			vector<vector<unsigned int> > const &dims,
			double const *lower, double const *upper,
			RNG *rng) const
{
    double sump = 0;
    unsigned int i = 0;

    for ( ; i < NCAT(dims); ++i) {
	sump += PROB(par)[i];
    }
    double p = sump * rng->uniform();    

    for ( ; i > 1; --i) {
	sump -= PROB(par)[i-1];
	if (sump <= p) 
            break;
    }
    *x  = i;
}

void DCat::support(double *lower, double *upper, unsigned int length,
	           vector<double const *> const &par,
	           vector<vector<unsigned int> > const &dims) const
{
    if (length != 1)
	throw logic_error("Invalid length in DCat::support");

    *lower = 1;
    *upper = NCAT(dims);
}

void DCat::typicalValue(double *x, unsigned int length,
			vector<double const *> const &par,
			vector<std::vector<unsigned int> > const &dims,
			double const *lower, double const *upper) const
{
    *x = max_element(PROB(par), PROB(par) + NCAT(dims)) - PROB(par) + 1;
}

bool DCat::isSupportFixed(std::vector<bool> const &fixmask) const
{
    return true;
}
