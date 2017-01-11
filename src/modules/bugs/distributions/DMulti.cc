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

bool DMulti::checkParameterLength(vector<unsigned int> const &len) const
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
			    vector<unsigned int> const &len) const
{
    if (SIZE(par) < 0)
	return false;

    // If SIZE is non-zero, we need at least one non-zero probability
    bool nz = SIZE(par) == 0;

    for (unsigned int i = 0; i < len[0]; ++i) {
	if (PROB(par)[i] < 0)
	    return false;
	else if (PROB(par)[i] > 0)
	    nz = true; 
    }

    return nz;
}

double DMulti::logDensity(double const *x, unsigned int length, PDFType type,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
    double loglik = 0.0;
    double S = 0;
    for (unsigned int i = 0; i < length; i++) {
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
	for (unsigned int i = 0; i < length; ++i) {
	    sump += PROB(par)[i];
	}
	loglik -= SIZE(par) * log(sump);
    }

    if (type != PDF_LIKELIHOOD) {
	//Terms depending on sampled value only
	for (unsigned int i = 0; i < length; ++i) {
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

void DMulti::randomSample(double *x, unsigned int length,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
    /* Sample multinomial as a series of binomial distributions */

    double N = SIZE(par);
    double const *prob = PROB(par);

    //Normalize probability
    double sump = 0;
    for (unsigned int i = 0; i < length; ++i) {
	sump += prob[i];
    }

    for (unsigned int i = 0; i < length - 1; i++) {
	if (N == 0) {
	    x[i] = 0;
	}
	else {
	    x[i] = rbinom(N, prob[i]/sump, rng);
	    N -= x[i];
	    sump -= prob[i];
	}
    }
    x[length - 1] = N;
}

void DMulti::support(double *lower, double *upper, unsigned int length,
	     vector<double const *> const &par,
	     vector<unsigned int> const &len) const
{
    for (unsigned int i = 0; i < length; ++i) {
	lower[i] = 0;
        if (PROB(par)[i] == 0) 
           upper[i] = 0;
        else
	   upper[i] = SIZE(par);
    }
}

unsigned int DMulti::length(vector<unsigned int> const &len) const
{
    return len[0];
}


void DMulti::typicalValue(double *x, unsigned int length,
			  vector<double const *> const &par,
			  vector<unsigned int> const &len,
			  double const *lower, double const *upper) const
{
    /* Draw a typical value in the same way as a random sample, but
       substituting the median at each stage */

    double N = SIZE(par);
    double const *prob = PROB(par);

    double sump = 0.0;
    for (unsigned int i = 0; i < length; ++i) {
	sump += prob[i];
    }

    for (unsigned int i = 0; i < length - 1; i++) {
	if (N == 0) {
	    x[i] = 0;
	}
	else {
	    x[i] = qbinom(0.5, N, prob[i]/sump, true, false);
	    N -= x[i];
	    sump -= prob[i];
	}
    }
    x[length - 1] = N;
}

bool DMulti::isSupportFixed(vector<bool> const &fixmask) const
{
    return fixmask[1];
}

unsigned int DMulti::df(vector<unsigned int> const &len) const
{
    return len[0] - 1;
} 

double DMulti::KL(vector<double const *> const &par1,
		  vector<double const *> const &par2,
		  vector<unsigned int> const &lengths) const
{
    if (SIZE(par1) != SIZE(par2))
	return JAGS_POSINF;

    unsigned int ncat = lengths[0];
    double y = 0, S1 = 0, S2 = 0;
    for (unsigned int i = 0; i < ncat; ++i) {
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
