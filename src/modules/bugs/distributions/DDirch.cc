#include <config.h>
#include <rng/RNG.h>
#include <sarray/SArray.h>
#include "DDirch.h"

#include <stdexcept>
#include <cmath>

#include <Rmath.h>

using std::vector;
using std::length_error;
using std::logic_error;

static inline double const *ALPHA(vector<SArray const *> const &par)
{
    return par[0]->value();
}

static inline unsigned long LENGTH(vector<SArray const *> const &par)
{
    return par[0]->length();
}

/* The Dirichlet distribution is extended to allow zero shape parameters.
   These represent structural zeros: when x ~ ddirch(alpha) is forward
   sampled, x[i] = 0 when alpha[i] = 0. To avoid trapping states in the
   model, structural zeros are only allowed when the alpha parameter
   is fixed.
*/

DDirch::DDirch()
  : Distribution("ddirch", 1, false, false) 
{}

DDirch::~DDirch()
{}

Index DDirch::dim(vector<Index> const &dims) const
{
    return dims[0];
}

bool DDirch::checkParameterDim(vector<Index> const &dims) const
{
    return (dims[0].size() == 1) && (dims[0][0] > 1);
}

bool DDirch::checkParameterValue(vector<SArray const *> const &par) const
{
    double const *alpha = ALPHA(par);
    unsigned long length = LENGTH(par);

    unsigned long nzero = 0; // No of zero shape par
    for (unsigned int i = 0; i < length; i++) {
	if (alpha[i] <= 0) {
	    if (alpha[i] == 0) {
		nzero++;
	    }
	    else {
		return false;
	    }
	}
    }

    if (nzero) {
	return par[0]->isFixed() && (nzero < length);
    }
    else {
	return true;
    }
}

double DDirch::logLikelihood(SArray const &x,
			     vector<SArray const *> const &par) const
{
    double const *y = x.value();
    double const *alpha = ALPHA(par);
    unsigned long length = LENGTH(par);

    double alphasum = 0.0;
    double loglik = 0.0;
    for (unsigned long i = 0; i < length; i++) {
	if (alpha[i] == 0) {
	    if (y[i] > 0)
		return -DBL_MAX;
	}
	else {
	    alphasum += alpha[i];
	    loglik += (alpha[i] - 1) * y[i] - lgamma(alpha[i]);
	}
    }
    loglik += lgamma(alphasum);

    return loglik;
}

void DDirch::randomSample(SArray &x, vector<SArray const *> const &par,
			  RNG *rng) const
{
    double const *alpha = ALPHA(par);
    unsigned long length = LENGTH(par);

    if (length != x.length()) 
	throw length_error ("Length mismatch error in DDirch::randomSample");

    /* Generate independent gamma random variables, then normalize
       to create Dirichlet distribution.
    */
    double *y = new double[length];
    double sumy = 0.0;
    for (unsigned long i = 0; i < length; i++) {
	y[i] = (alpha[i]==0) ? 0 : rgamma(alpha[i], 1, rng);
	sumy += y[i];
    }
    for (unsigned long j = 0; j < length; j++) {
	y[j] /= sumy;
    }
    x.setValue(y, length);
    delete [] y;
}

unsigned long DDirch::df(std::vector<SArray const *> const &par) const
{
    unsigned long length = LENGTH(par);
    double const *alpha = ALPHA(par);

    unsigned long d = length - 1;
    for (unsigned long i = 0; i < length; ++i) {
	if (d == 0) {
	    throw logic_error("Bad degrees of freedom in DDirch");
	}
	else {
	    d -= (alpha[i] == 0);
	}
    }
    return d;
}

double 
DDirch::lowerSupport(unsigned long i,
		     std::vector<SArray const *> const &par) const
{
    if (i >= LENGTH(par))
	throw logic_error("Invalid index in DDirch::lowerSupport");
    
    return 0;
}

double 
DDirch:: upperSupport(unsigned long i,
		      std::vector<SArray const *> const &par) const
{
    if (i >= LENGTH(par))
	throw logic_error("Invalid index in DDirch::upperSupport");

    if (ALPHA(par)[i] == 0)
	return 0;
    else
	return DBL_MAX;

}

void DDirch::typicalValue(SArray &x, vector<SArray const *> const &par) const
{
  unsigned long N = LENGTH(par);

  double alphasum = 0.0;
  for (unsigned long i = 0; i < N; ++i) {
    alphasum += ALPHA(par)[i];
  }
  for (unsigned long i = 0; i < N; ++i) {
    x.setValue(ALPHA(par)[i]/alphasum, i);
  }
}
