#include <config.h>

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <module/ModuleError.h>
#include <util/integer.h>

#include "DScaledWishart.h"

#include <cfloat>
#include <cmath>
#include <vector>
#include <algorithm>

#include <JRmath.h>

using std::vector;
using std::log;
using std::copy;

#define SCALE(par) (par[0])
#define DF(par)    (*par[1])
#define NROW(dims)  (dims[0][0])

#define F77_DSYEV F77_FUNC(dsyev,DSYEV)

extern "C" {
    void F77_DSYEV (const char* jobz, const char* uplo,
		    const int* n, double* a, const int* lda,
		    double* w, 
		    double* work, const int* lwork, int* info);
}

namespace jags {
namespace glm {

    //FIXME: This is copy-pasted from the bugs module
    static double logdet(double const *a, unsigned long n)
    {
	// Log determinant of n x n symmetric positive matrix a */
  
	unsigned long N = n*n;
	vector<double> acopy(N);
	copy(a, a+N, acopy.begin());

	vector<double> w(n);
	int lwork = -1;
	double worktest = 0;
	int info = 0;
	int ni = asInteger(n);
	F77_DSYEV("N","L", &ni, &acopy[0], &ni, &w[0], &worktest, &lwork, &info);
	if (info != 0) {
	    throwRuntimeError("unable to calculate workspace size for dsyev");
	}
	lwork = static_cast<int>(worktest);
	double *work = new double[lwork];
	F77_DSYEV("N","L", &ni, &acopy[0], &ni, &w[0], work, &lwork, &info);
	delete [] work;
	if (info != 0) {
	    throwRuntimeError("unable to calculate eigenvalues in dsyev");
	}

	if (w[0] <= 0) {
	    throwRuntimeError("Non positive definite matrix in call to logdet");
	}

	double logdet = 0;
	for (unsigned long i = 0; i < n; i++) {
	    logdet += log(w[i]);
	}
	
	return logdet;
    }
    
    static double log_multigamma(double n, unsigned long p)
    {
	double y =  (p * (p-1) * log(M_PI))/4;
	for (unsigned long j = 0; j < p; ++j) {
	    y += lgammafn((n-j)/2);
	}
	return y;
    }
    
    DScaledWishart::DScaledWishart()
	: ArrayDist("dscaled.wishart", 2) 
    {}

    double
    DScaledWishart::logDensity(double const *x, PDFType type,
			       vector<double const *> const &par,
			       vector<vector<unsigned long> > const &dims) const
    {
	double const *A = SCALE(par);
	unsigned long p = NROW(dims);
	double df = DF(par);
	double k = p + df - 1;
	
	double loglik = (k - p - 1) * logdet(x, p) / 2;
	for (unsigned int i = 0; i < p; ++i) {
	    loglik -= (k+1) * log(df * x[i*p + i] + 1/(A[i]*A[i])) / 2;
	    
	}

	if (type != PDF_PRIOR) {
	    //Normalize density
	    loglik += p * k * log(df) / 2;
	    for (unsigned int i = 0; i < p; ++i) {
		loglik -= log(A[i]);
	    }
	    loglik += p * lgammafn((k+1)/2) -  p * lgammafn(0.5) -
		- log_multigamma(k, p);
	}

	return loglik;
    }

    void DScaledWishart::sampleWishart(double *x,
				       double const *R, unsigned long nrow,
				       double k, RNG *rng)
    {
	/* 
	   Generate random Wishart variable with diagonal prior matrix
	   Equivalent to DWish(diag(R), k)
	   
	   The algorithm is adapted from bugs::DWish::randomSample but
	   is considerably simpler.
	*/

	/* Generate square root of Wishart random variable:
	   - diagonal elements are square root of Chi square
	   - upper off-diagonal elements are normal
	   - lower off-diagonal elements are zero
	*/
	vector<vector<double> > Z(nrow, vector<double>(nrow));

	for (unsigned int i = 0; i < nrow; ++i) {
	    for (unsigned int l = 0; l < i; ++l) {
		Z[i][l] = rnorm(0, 1, rng);
	    }
	    Z[i][i] = sqrt(rchisq(k - i, rng));    
	}

	/* Take inverse square root of R */
	vector<double> scale(nrow);
	for (unsigned int i = 0; i < nrow; ++i) {
	    scale[i] = 1/sqrt(R[i]);
	}
	
	/* Now put cross-product into x */
	for (unsigned int i = 0; i < nrow; i++) {
	    for (unsigned int j = 0; j <= i; j++) {
		double xx = 0;
		for (unsigned int l = 0; l <= j; l++) {
		    xx += Z[i][l] * Z[j][l];
		}
		x[nrow * j + i] = x[nrow * i + j] = scale[j] * scale[i] * xx;
	    }
	}
    }

void DScaledWishart::randomSample(double *x,
				  vector<double const *> const &par,
				  vector<vector<unsigned long> > const &dims,
				  RNG *rng) const
{
    unsigned long nrow = NROW(dims);
    double df = DF(par);
    double const *scale = SCALE(par);
    
    vector<double> R(nrow);
    double k = nrow + df - 1;
    for (unsigned int i = 0; i < nrow; ++i) {
	//NB Rmath implementation of gamma is parameterized in terms of
	//shape and scale, unlike BUGS which uses shape and rate
	R[i] = 2 * df * rgamma(0.5, scale[i]*scale[i], rng);
    }
    sampleWishart(x, &R[0], nrow, k, rng);
}

bool
DScaledWishart::checkParameterDim (vector<vector<unsigned long> > const &dims)
    const
{
    return (isVector(dims[0]) || isScalar(dims[0])) && isScalar(dims[1]);
}

vector<unsigned long> 
DScaledWishart::dim(vector<vector<unsigned long> > const &dims) const
{
    if (isScalar(dims[0])) {
	return vector<unsigned long>(1,1);
    }
    else {
	return vector<unsigned long> (2, dims[0][0]);
    }
}

bool 
DScaledWishart::checkParameterValue(vector<double const *> const &par,
			   vector<vector<unsigned long> > const &dims) const
{
    if (DF(par) < 1) return false;
    double const *scale = SCALE(par);
    unsigned long n = NROW(dims);
    for (unsigned long i = 0; i < n; ++i) {
	if (scale[i] <= 0) return false;
    }
    return true;
}


void DScaledWishart::support(double *lower, double *upper,
		    vector<double const *> const &par,
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
    
bool DScaledWishart::isSupportFixed(vector<bool> const &fixmask) const
{
    return true;
}

    bool DScaledWishart::fullRank() const
    {   
	return false;
    }

}}
