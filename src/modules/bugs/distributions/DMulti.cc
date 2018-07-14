#include <config.h>
#include "DMulti.h"
#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>

#include <cmath>

#include <JRmath.h>

using std::vector;
using std::string;

#define PROB(par) (par[0])
#define SIZE(par) (*par[1])

namespace jags {
namespace bugs {

DMulti::DMulti()
  : VectorDist("dmulti", 2) 
{}

string DMulti::alias() const
{
    return "dmultinom";
}
    
bool DMulti::isDiscreteValued(vector<bool> const &mask) const
{
    return true;
}

bool DMulti::checkParameterLength(vector<unsigned long> const &len) const
{
    //Check that PROB is non-empty and SIZE is a scalar
    return len[0] >= 1 && len[1] == 1;
}

bool DMulti::checkParameterDiscrete(vector<bool> const &mask) const
{
    return mask[1]; //SIZE is discrete-valued
}

bool 
DMulti::checkParameterValue(vector<double const *> const &par,
			    vector<unsigned long> const &len) const
{
    if (SIZE(par) < 0)
	return false;

    // If SIZE is non-zero, we need at least one non-zero probability
    bool nz = SIZE(par) == 0;

    for (unsigned long i = 0; i < len[0]; ++i) {
	if (PROB(par)[i] < 0)
	    return false;
	else if (PROB(par)[i] > 0)
	    nz = true; 
    }

    return nz;
}

double DMulti::logDensity(double const *x, PDFType type,
			  vector<double const *> const &par,
			  vector<unsigned long> const &len,
			  double const *lower, double const *upper) const
{
    double loglik = 0.0;
    double S = 0;
    unsigned long N = length(len);
    for (unsigned long i = 0; i < N; i++) {
	if (x[i] < 0 || floor(x[i]) != x[i]) {
	    return JAGS_NEGINF;
	}
	else if (x[i] != 0) {
	    if (PROB(par)[i] == 0) {
		return JAGS_NEGINF;
	    }
	    else {
		loglik += x[i] * log(PROB(par)[i]);
	    }
	    S += x[i];
	}
    }

    //Check consistency between parameters and data
    if (S != SIZE(par))
	return JAGS_NEGINF;

    if (type != PDF_PRIOR) {
	//Terms depending on parameters only
	double sump = 0.0;
	for (unsigned long i = 0; i < N; ++i) {
	    sump += PROB(par)[i];
	}
	if (SIZE(par) != 0) {
	    loglik -= SIZE(par) * log(sump);
	}
    }

    if (type != PDF_LIKELIHOOD) {
	//Terms depending on sampled value only
	for (unsigned long i = 0; i < N; ++i) {
	    loglik -= lgammafn(x[i] + 1);
	}
    }

    if (type == PDF_FULL) {
	//If either data or parameters are fixed then this term is constant
	//bearing in mind consistency check above.
	loglik += lgammafn(SIZE(par) + 1);
    }

    return loglik;
}

void DMulti::randomSample(double *x,
			  vector<double const *> const &par,
			  vector<unsigned long> const &len,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
    /* Sample multinomial as a series of binomial distributions */

    double size = SIZE(par);
    double const *prob = PROB(par);
    unsigned long N = length(len);
    
    //Normalize probability
    double sump = 0;
    for (unsigned long i = 0; i < N; ++i) {
	sump += prob[i];
    }

    for (unsigned long i = 0; i < N - 1; i++) {
	if (size == 0) {
	    x[i] = 0;
	}
	else {
	    x[i] = rbinom(size, prob[i]/sump, rng);
	    size -= x[i];
	    sump -= prob[i];
	}
    }
    x[N - 1] = size;
}

void DMulti::support(double *lower, double *upper,
	     vector<double const *> const &par,
	     vector<unsigned long> const &len) const
{
    unsigned long N = length(len);
    for (unsigned long i = 0; i < N; ++i) {
	lower[i] = 0;
        if (PROB(par)[i] == 0) 
           upper[i] = 0;
        else
	   upper[i] = SIZE(par);
    }
}

unsigned long DMulti::length(vector<unsigned long> const &len) const
{
    return len[0];
}

bool DMulti::isSupportFixed(vector<bool> const &fixmask) const
{
    return fixmask[1];
}

unsigned long DMulti::df(vector<unsigned long> const &len) const
{
    return len[0] - 1;
} 

double DMulti::KL(vector<double const *> const &par1,
		  vector<double const *> const &par2,
		  vector<unsigned long> const &lengths) const
{
    if (SIZE(par1) != SIZE(par2))
	return JAGS_POSINF;

    unsigned long ncat = lengths[0];
    double y = 0, S1 = 0, S2 = 0;
    for (unsigned long i = 0; i < ncat; ++i) {
	double p1 = PROB(par1)[i];
	double p2 = PROB(par2)[i];
	
	if (p1 == 0) {
	    S2 += p2;
	}
	else if (p2 == 0) {
	    return JAGS_POSINF;
	}
	else {
	    y += p1 * (log(p1) - log(p2));
	    S1 += p1;
	    S2 += p2;
	}
    }
    y /= S1;
    y += log(S2) - log(S1);
    y *= SIZE(par1);

    return y;
}


}}
