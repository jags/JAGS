#include <config.h>
#include "DRW1.h"
#include <util/dim.h>
#include <util/nainf.h>

#include <cmath>

#include <JRmath.h>
#include <rng/RNG.h>

using std::vector;
using std::string;
using std::abs;

#define TAU(par) (*par[0])
#define XCOORDS(par) (par[1])

namespace jags {
    namespace bugs {

	DRW1::DRW1()
	    : VectorDist("drw1", 2) 
	{}

	unsigned long DRW1::length(vector<unsigned long> const &len) const
	{
	    return len[1];
	}

	bool DRW1::checkParameterLength(vector<unsigned long> const &len) const
	{
	    return len[0] == 1 && len[1] > 1;
	}

	bool
	DRW1::checkParameterValue(vector<double const *> const &par,
				  vector<unsigned long> const &len) const
	{
	    unsigned long N = len[1];
	    double const *xcoords = XCOORDS(par);
	    double tau = TAU(par);

	    //Only regular random walks are supported
	    for (unsigned long i = 1; i < N; ++i) {
		double delta = xcoords[i] - xcoords[i-1];
		if (abs(delta - 1) > 1.0E-6) return false;
	    }

	    if (tau < 0) return false;

	    return true;
	}

	double DRW1::logDensity(double const *x, unsigned long length,
				PDFType ,
				vector<double const *> const &par,
				vector<unsigned long> const &len,
				double const *, double const *) const
	{
	    double tau = TAU(par);
	    double sigma = 1/sqrt(tau);
	    
	    double loglik = 0.0;
	    for (unsigned long i = 1; i < length; i++) {
		loglik += dnorm(x[i] - x[i-1], 0, sigma, true);
	    }
	    return loglik;
	}

	void DRW1::randomSample(double *x, unsigned long length,
				vector<double const *> const &par,
				vector<unsigned long> const &len,
				double const *, double const *,
				RNG *rng) const
	{
	    double tau = TAU(par);
	    double sigma = 1/sqrt(tau);

	    double S = 0.0;
	    x[0] = 0;
	    for (unsigned long i = 1; i < length; i++) {
		x[i] = x[i-1] + rnorm(0, sigma, rng);
		S += x[i];
	    }
	    S /= length;
	    for (unsigned long i = 0; i < length; i++) {
		x[i] -= S;
	    }
	}

	void DRW1::support(double *lower, double *upper, unsigned long length,
			   vector<double const *> const &,
			   vector<unsigned long> const &len) const
	{
	    for (unsigned long i = 0; i < length; ++i) {
		lower[i] = JAGS_NEGINF;
		upper[i] = JAGS_POSINF;
	    }
	}

	bool DRW1::isSupportFixed(vector<bool> const &) const
	{
	    return true;
	}

	unsigned long DRW1::df(vector<unsigned long> const &len) const
	{
	    return len[1] - 1;
	}

    }
}
