#include <config.h>

#include <rng/RNG.h>
#include <sarray/SArray.h>
#include <sarray/util.h>

#include "lapack.h"
#include "matrix.h"
#include "DWish.h"

#include <stdexcept>
#include <cfloat>
#include <cmath>
#include <vector>

#include <Rmath.h>

using std::vector;
using std::logic_error;
using std::runtime_error;

static inline double const *SCALE(vector<SArray const *> const &par)
{
    return par[0]->value();
}

static inline double DF(vector<SArray const *> const &par)
{
    return  *par[1]->value();
}

static inline unsigned long NROW(vector<SArray const *> const &par)
{
    return par[0]->dim(true)[0];
}

static inline unsigned long LENGTH(vector<SArray const *> const &par)
{
    return par[0]->length();
}

DWish::DWish()
  : Distribution("dwish", 2, false, false) 
{}

DWish::~DWish()
{}

double DWish::logLikelihood(SArray const &x,
			    vector<SArray const *> const &par) const
{
    double const *y = x.value();
    double const *scale = SCALE(par);
    long p = NROW(par);

    double loglik = 0;
    long length = LENGTH(par);
    for (int i = 0; i < length; ++i) {
	loglik += scale[i] * y[i];
    }
    loglik += DF(par) * logdet(scale, p) + (DF(par) - p - 1) * logdet(y, p);
    return loglik/2;
}


void DWish::randomSample(SArray &x, vector<SArray const *> const &par,
			 RNG *rng) const
{
  long length = LENGTH(par);
  double *y = new double[length];
  randomSample(y, length, SCALE(par), DF(par), NROW(par), rng);
  x.setValue(y, length);
  delete [] y;
}

void DWish::randomSample(double *x, long length,
			 double const *R, double k, int nrow,
                         RNG *rng)
{
    /* 
       Generate random Wishart variable, using an algorithm proposed
       by Bill Venables and originally implemented in S.
    */
    if (length != nrow*nrow) {
	throw logic_error("invalid length in DWish::randomSample");
    }

    /* 
       Get inverse of R. Venables' algorithm was implemented in
       terms of the inverse of R, but we use a different parameterization
       to preserve conjugacy.
    */
    double * C = new double[length];
    inverse(C, R, nrow, true);

    /* Get Choleskly decomposition of C */
    int info = 0;
    F77_DPOTRF("U", &nrow, C, &nrow, &info);
    if (info != 0) {
	throw runtime_error("Failed to get Cholesky decomposition of R in dwish");
    }
    
    /* Set lower triangle of C to zero */
    for (int j = 0; j < nrow; j++) {
	double * C_j = &C[j*nrow]; //column j of matrix C
	for (int i = j + 1; i < nrow; i++) {
	    C_j[i] = 0;
	}
    }

    /* Generate square root of Wishart random variable:
       - diagonal elements are square root of Chi square
       - upper off-diagonal elements are normal
       - lower off-diagonal elements are zero
    */
    double *Z = new double[length];
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
  
    /* Transform Z with Cholesky decomposition */
    double *Ztrans = new double[length];
    for (int i = 0; i < nrow; i++) {
	for (int j = 0; j < nrow; j++) {
	    double zz = 0;
	    for (int l = 0; l < nrow; l++) {
		zz += Z[nrow * l + i] * C[nrow * j + l];
	    }
	    Ztrans[nrow * j + i] = zz;
	}
    }
    delete [] C;
    delete [] Z;

    /* Now put cross-product into x */
    for (int i = 0; i < nrow; i++) {
	double const *Ztrans_i = &Ztrans[nrow * i];
	for (int j = 0; j <= i; j++) {
	    double const *Ztrans_j = &Ztrans[nrow * j];
	    double xx = 0;
	    for (int l = 0; l < nrow; l++) {
		xx += Ztrans_i[l] * Ztrans_j[l];
	    }
	    x[nrow * j + i] = x[nrow * i + j] = xx;
	}
    }
    delete [] Ztrans;
}

bool DWish::checkParameterDim (vector<Index> const &dims) const
{
  return isSquareMatrix(dims[0]) && isScalar(dims[1]);
}

Index DWish::dim(vector<Index> const &dims) const
{
  return dims[0];
}

bool DWish::checkParameterValue (vector<SArray const *> const &par) const
{
    /* Check that we have sufficient degrees of freedom */
    if (DF(par) < NROW(par))
	return false;
    /* Check symmetry of scale matrix */
    double const *scale = SCALE(par);
    long nrow = NROW(par);
    for (int i = 0; i < nrow; ++i) {
	for (int j = 0; j < i; ++i) {
	    if (fabs(scale[i + nrow*j] - scale[j + nrow*i]) > DBL_EPSILON)
		return false;
	}
    }
    /* Skipping check of positive definiteness of scale matrix */
    return true;
}

unsigned long DWish::df(std::vector<SArray const *> const &par) const
{
    return LENGTH(par);
}

double 
DWish::lowerSupport(unsigned long i,
		    std::vector<SArray const *> const &par) const
{
    if (i % NROW(par) == i / NROW(par)) {
	//Diagonal elements
	return 0;
    }
    else {
	return -DBL_MAX;
    }
}

double 
DWish::upperSupport(unsigned long i,
		    std::vector<SArray const *> const &par) const
{
    return DBL_MAX;
}

void 
DWish::typicalValue(SArray &x, std::vector<SArray const *> const &par)
  const
{
  /* Returns the mean as a typical value. We need to invert the
     scale matrix */

  unsigned long length = LENGTH(par);
  double * C = new double[length];
  inverse(C, SCALE(par), NROW(par), true);
  for (unsigned long i = 0; i < length; ++i) {
    C[i] *= DF(par);
  }
  x.setValue(C, length);
  delete [] C;
}
