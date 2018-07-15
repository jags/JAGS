#include <config.h>

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <module/ModuleError.h>
#include <util/integer.h>

#include "lapack.h"
#include "matrix.h"
#include "DWish.h"

#include <cfloat>
#include <cmath>
#include <vector>
#include <algorithm>

#include <JRmath.h>

using std::vector;
using std::log;
using std::reverse;

#define SCALE(par) (par[0])
#define DF(par)    (*par[1])
#define NROW(dims)  (dims[0][0])

static double log_multigamma(double n, unsigned long p)
{
    double y =  (p * (p-1) * log(M_PI))/4;
    for (unsigned long j = 0; j < p; ++j) {
        y += lgammafn((n-j)/2);
    }
    return y;
}

namespace jags {
namespace bugs {
    
DWish::DWish()
  : ArrayDist("dwish", 2) 
{}

double DWish::logDensity(double const *x, PDFType type,
			 vector<double const *> const &par,
			 vector<vector<unsigned long> > const &dims) const
{
    double const *scale = SCALE(par);
    unsigned long p = NROW(dims);
    unsigned long length = p * p;
    
    double loglik = (DF(par) - p - 1) * logdet(x, p);
    for (unsigned long i = 0; i < length; ++i) {
	loglik -= scale[i] * x[i];
    }

    if (type != PDF_PRIOR) {
	//Normalize density
	loglik += DF(par) * logdet(scale, p) -  DF(par) * p * log(2.0) -
	    2 * log_multigamma(DF(par), p);
    }

    return loglik/2;
}

void DWish::randomSample(double *X,
			 double const *R, double k, unsigned long nrow,
			 RNG *rng)
{
    /* 
       Generate random Wishart variable, using an algorithm proposed
       by Bill Venables and originally implemented in S.
    */
    unsigned long length = nrow * nrow;

    /* 
       Get Cholesky decomposition of the inverse of R. First we
       factorize R (dpotrf) and then we invert the triangular factor
       (dtrtri). At the end, C contains the upper triangular Cholesky
       factor. NB We must reverse the elements of the matrix at start
       and end of the calculations.
    */
    vector<double> C(length);
    copy(R, R + length, C.rbegin());
    int info = 0;
    int ni = asInteger(nrow);
    F77_DPOTRF("L", &ni, &C[0], &ni, &info);
    if (info != 0) {
	jags::throwRuntimeError("Failed to get Cholesky decomposition of R");
    }
    F77_DTRTRI("L", "N", &ni, &C[0], &ni, &info);
    if (info != 0) {
	jags::throwRuntimeError("Failed to invert Cholesky decomposition of R");
    }
    reverse(C.begin(), C.end());
    
    /* Generate square root of Wishart random variable:
       - diagonal elements are square root of Chi square
       - upper off-diagonal elements are normal
       - lower off-diagonal elements are zero
    */
    vector<double> Z(length);
    for (unsigned long j = 0; j < nrow; j++) {
	double *Z_j = &Z[j*nrow]; //jth column of Z
	for (unsigned long i = 0; i < j; i++) {
	    Z_j[i] = rnorm(0, 1, rng);
	}
	Z_j[j] = sqrt(rchisq(k - j, rng));    
	for (unsigned long i = j + 1; i < nrow; i++) {
	    Z_j[i] = 0;
	}
    }

    // Z = Z %*% C 
    double one = 1;
    F77_DTRMM("R", "U", "N", "N", &ni, &ni, &one, &C[0], &ni, &Z[0], &ni);

    // X = t(Z) %*% Z
    double zero = 0;
    F77_DSYRK("U", "T", &ni, &ni, &one, &Z[0], &ni, &zero, X, &ni);

    // Copy lower triangle of X from upper triangle
    for (unsigned long i = 0; i < nrow; ++i) {
	for (unsigned long j = 0; j < i; ++j) {
	    X[j * nrow + i] = X[i * nrow + j];
	}
    }
}

void DWish::randomSample(double *x,
			 vector<double const *> const &par,
			 vector<vector<unsigned long> > const &dims,
			 RNG *rng) const
{
    randomSample(x, SCALE(par), DF(par), NROW(dims), rng);
}

bool DWish::checkParameterDim (vector<vector<unsigned long> > const &dims) const
{
  return isSquareMatrix(dims[0]) && isScalar(dims[1]);
}

vector<unsigned long> 
DWish::dim(vector<vector<unsigned long> > const &dims) const
{
  return dims[0];
}

bool 
DWish::checkParameterValue(vector<double const *> const &par,
			   vector<vector<unsigned long> > const &dims) const
{
    // Check that we have sufficient degrees of freedom
    return DF(par) >= NROW(dims);
}


void DWish::support(double *lower, double *upper,
		    vector<double const *> const &,
		    vector<vector<unsigned long> > const &dims) const
{
    unsigned long length = dims[0][0] * dims[0][1];
    for (unsigned long i = 0; i < length; ++i) {
	if (i % NROW(dims) == i / NROW(dims)) {
	    //Diagonal elements
	    lower[i] =  0;
	}
	else {
	    lower[i] =  JAGS_NEGINF;
	}
	upper[i] = JAGS_POSINF;
    }
}

bool DWish::isSupportFixed(vector<bool> const &) const
{
    return true;
}

unsigned long DWish::df(vector<vector<unsigned long> > const &dims) const
{   
  return dims[0][0] * (dims[0][0] + 1) / 2;
}

}}
