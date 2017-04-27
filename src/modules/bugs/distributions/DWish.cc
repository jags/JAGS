#include <config.h>

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <module/ModuleError.h>

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

static double log_multigamma(double n, unsigned int p)
{
    double y =  (p * (p-1) * log(M_PI))/4;
    for (unsigned int j = 0; j < p; ++j) {
        y += lgammafn((n-j)/2);
    }
    return y;
}

namespace jags {
namespace bugs {
    
DWish::DWish()
  : ArrayDist("dwish", 2) 
{}

double DWish::logDensity(double const *x, unsigned int length, PDFType type,
			 vector<double const *> const &par,
			 vector<vector<unsigned int> > const &dims,
			 double const *lower, double const *upper) const
{
    double const *scale = SCALE(par);
    unsigned int p = NROW(dims);

    double loglik = (DF(par) - p - 1) * logdet(x, p);
    for (unsigned int i = 0; i < length; ++i) {
	loglik -= scale[i] * x[i];
    }

    if (type != PDF_PRIOR) {
	//Normalize density
	loglik += DF(par) * logdet(scale, p) -  DF(par) * p * log(2.0) -
	    2 * log_multigamma(DF(par), p);
    }

    return loglik/2;
}

void DWish::randomSample(double *X, int length,
			 double const *R, double k, int nrow,
			 RNG *rng)
{
    /* 
       Generate random Wishart variable, using an algorithm proposed
       by Bill Venables and originally implemented in S.
    */

    if (length != nrow*nrow) {
	jags::throwLogicError("invalid length in DWish::randomSample");
    }

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
    F77_DPOTRF("L", &nrow, &C[0], &nrow, &info);
    if (info != 0) {
	jags::throwRuntimeError("Failed to get Cholesky decomposition of R");
    }
    F77_DTRTRI("L", "N", &nrow, &C[0], &nrow, &info);
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
    for (int j = 0; j < nrow; j++) {
	double *Z_j = &Z[j*nrow]; //jth column of Z
	for (int i = 0; i < j; i++) {
	    Z_j[i] = rnorm(0, 1, rng);
	}
	Z_j[j] = sqrt(rchisq(k - j, rng));    
	for (int i = j + 1; i < nrow; i++) {
	    Z_j[i] = 0;
	}
    }

    // Z = Z %*% C 
    double one = 1;
    F77_DTRMM("R", "U", "N", "N", &nrow, &nrow, &one, &C[0], &nrow, &Z[0],
	      &nrow);

    // X = t(Z) %*% Z
    double zero = 0;
    F77_DSYRK("U", "T", &nrow, &nrow, &one, &Z[0], &nrow, &zero, X, &nrow);

    // Copy lower triangle of X from upper triangle
    for (int i = 0; i < nrow; ++i) {
	for (int j = 0; j < i; ++j) {
	    X[j * nrow + i] = X[i * nrow + j];
	}
    }
}

void DWish::randomSample(double *x, unsigned int length,
			 vector<double const *> const &par,
			 vector<vector<unsigned int> > const &dims,
			 double const *lower, double const *upper,
			 RNG *rng) const
{
    randomSample(x, length, SCALE(par), DF(par), NROW(dims), rng);
}

bool DWish::checkParameterDim (vector<vector<unsigned int> > const &dims) const
{
  return isSquareMatrix(dims[0]) && isScalar(dims[1]);
}

vector<unsigned int> 
DWish::dim(vector<vector<unsigned int> > const &dims) const
{
  return dims[0];
}

bool 
DWish::checkParameterValue(vector<double const *> const &par,
			   vector<vector<unsigned int> > const &dims) const
{
    // Check that we have sufficient degrees of freedom
    if (DF(par) < NROW(dims)) return false;
    // Check scale matrix

    double const *scale = SCALE(par);
    unsigned int n = NROW(dims);
    return check_symmetry(scale, n) && check_symmetric_ispd(scale, n);
}


void DWish::support(double *lower, double *upper, unsigned int length,
		    vector<double const *> const &par,
		    vector<vector<unsigned int> > const &dims) const
{
    for (unsigned int i = 0; i < length; ++i) {
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

void DWish::typicalValue(double *x, unsigned int length,
			 vector<double const *> const &par,
			 vector<vector<unsigned int> > const &dims,
			 double const *lower, double const *upper) const
{
    /* Returns the mean as a typical value. We need to invert the
       scale matrix */

    if (!inverse_spd(x, SCALE(par), NROW(dims))) {
	throwDistError(this, "Inverse failed in typicalValue");
    }
    for (unsigned int i = 0; i < length; ++i) {
	x[i] *= DF(par);
    }
}

bool DWish::isSupportFixed(vector<bool> const &fixmask) const
{
    return true;
}

unsigned int DWish::df(vector<vector<unsigned int> > const &dims) const
{   
  return dims[0][0] * (dims[0][0] + 1) / 2;
}

}}
