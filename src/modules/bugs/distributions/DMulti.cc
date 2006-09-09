#include <config.h>
#include <rng/RNG.h>
#include <sarray/util.h>
#include <sarray/SArray.h>
#include "DMulti.h"

#include <stdexcept>
#include <cmath>

#include <Rmath.h>

using std::vector;
using std::length_error;
using std::logic_error;

static inline double const *PROB(vector<SArray const *> const &parameters)
{
    return parameters[0]->value();
}

static inline unsigned long LENGTH(vector<SArray const *> const &parameters)
{
    return parameters[0]->length();
}

static inline double SIZE(vector<SArray const *> const &parameters) 
{
    return *parameters[1]->value();
}

DMulti::DMulti()
  : Distribution("dmulti", 2, false, true) 
{}

DMulti::~DMulti()
{}

bool DMulti::checkParameterDim(vector<Index> const &dims) const
{
  //Check that PROB is a vector and SIZE is a scalar
  return isVector(dims[0]) && isScalar(dims[1]);
}

bool DMulti::checkParameterDiscrete(vector<bool> const &mask) const
{
    return mask[1]; //SIZE is discrete-valued
}

bool DMulti::checkParameterValue(vector<SArray const *> const &par) const
{
    double const *prob = PROB(par);
    unsigned long length = LENGTH(par);
  
    if (SIZE(par) < 1)
	return false;

    double sump = 0.0;
    for (unsigned long i = 0; i < length; i++) {
	if (prob[i] < 0 || prob[i] > 1)
	    return false;
	sump += prob[i];
    }
    if (fabs(sump - 1.0) > 1.0E-6)
	return false;

    return true;
}

double DMulti::logLikelihood(SArray const &x,
			     vector<SArray const *> const &par) const
{
    double const *y = x.value();
    double const *prob = PROB(par);
    unsigned long length = LENGTH(par);

    double loglik = 0.0;
    double ysum = 0.0;
    for (unsigned long i = 0; i < length; i++) {
	if (prob[i] == 0) {
	    if (y[i] != 0)
		return -DBL_MAX;
	}
	else if (y[i] != 0) {
	     loglik += y[i] * log(prob[i]) - lgamma(y[i] + 1);
	     ysum += y[i];
	}
    }
    loglik += lgamma(ysum + 1);

    return loglik;
}

void DMulti::randomSample(SArray &x, vector<SArray const *> const &par,
			  RNG *rng) const
{
    unsigned long length = LENGTH(par);

    /* Sample multinomial as a series of binomial distributions */
    double *y = new double[length];
    double N = SIZE(par);
    double sump = 1.0;
    double const *prob = PROB(par);
    for (unsigned long i = 0; i < length - 1; i++) {
	if (N == 0) {
	    y[i] = 0;
	}
	else {
	    y[i] = rbinom(N, prob[i]/sump, rng);
	    N -= y[i];
	    sump -= prob[i];
	}
    }
    y[length - 1] = N;
    x.setValue(y, length);
    delete [] y;
}

unsigned long DMulti::df(std::vector<SArray const *> const &par) const
{
    unsigned long length = LENGTH(par);
    double const *prob = PROB(par);
    //FIXME: Should we allow structural zeros only for fixed p, as in DDirch?
    unsigned long df = length - 1;
    for (unsigned long i = 0; i < length; ++i) {
	if (df == 0) {
	    throw logic_error("Bad degrees of freedom in DMulti");
	}
	else {
	    df -= (prob[i] == 0);
	}
    }
    return df;
}

double 
DMulti::lowerSupport(unsigned long i,
		     std::vector<SArray const *> const &par) const
{
    if (i >= LENGTH(par))
	throw logic_error("Invalid index in DMulti::lowerSupport");
    
    return 0.0;
}

double 
DMulti::upperSupport(unsigned long i,
		     std::vector<SArray const *> const &par) const
{
    if (i >= LENGTH(par))
	throw logic_error("Invalid index in DMulti::upperSupport");
    
    return DBL_MAX;
}

Index DMulti::dim(std::vector<Index> const &dims) const
{
    return dims[0];
}


void DMulti::typicalValue(SArray &x, vector<SArray const *> const &par) const
{
  /* Draw a typical value in the same way as a random sample, but
     substituting the median at each stage */

  unsigned long length = LENGTH(par);

  double *y = new double[length];
  double N = SIZE(par);
  double sump = 1.0;
  double const *prob = PROB(par);
  for (unsigned long i = 0; i < length - 1; i++) {
    if (N == 0) {
      y[i] = 0;
    }
    else {
      y[i] = qbinom(0.5, N, prob[i]/sump, true, false);
      N -= y[i];
      sump -= prob[i];
    }
  }
  y[length - 1] = N;
  x.setValue(y, length);
  delete [] y;
}
