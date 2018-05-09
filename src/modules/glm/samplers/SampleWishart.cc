#include <config.h>

#include "SampleWishart.h"

#include <rng/RNG.h>
#include <module/ModuleError.h>
#include <util/integer.h>

#include <vector>
#include <cmath>
#include <algorithm>

#include <JRmath.h>

using std::vector;
using std::sqrt;
using std::copy;
using std::reverse;


#define F77_DPOTRF F77_FUNC(dpotrf,DPOTRF)
#define F77_DTRTRI F77_FUNC(dtrtri, DTRTRI)
#define F77_DTRMM  F77_FUNC(dtrmm, DTRMM)
#define F77_DSYRK  F77_FUNC(dsyrk, DSYRK)

extern "C" {
    void F77_DPOTRF (const char *uplo, const int *n, double *a,
		     const int *lda, const int *info);
    void F77_DTRTRI (const char *uplo, const char *diag,
		     const int *n, double *a, const int *lda, const int *info);
    void F77_DTRMM(const char *side, const char *uplo, const char *transa,
		   const char *diag, const int *m, const int *n,
		   const double *alpha, const double *a, const int *lda,
		   double *b, const int *ldb);
    void F77_DSYRK(const char *uplo, const char *trans, const int *n,
		   const int *k,
		   const double *alpha, const double *a, const int *lda,
		   const double *beta, double *c, const int *ldc);
}


namespace jags {
    namespace glm {

	//FIXME We would not need this if we could call
	//bugs::DWish::sampleWishart
	void sampleWishart(double *X, unsigned long length,
			   double const *R, double df, unsigned long nrow,
			   RNG *rng)
	{
	    if (df <= nrow) {
		throwLogicError("Invalid df in sampleWishart");
	    }
	    
	    int info = 0;
	    /* 
	       Generate random Wishart variable, using an algorithm
	       proposed by Bill Venables and originally implemented in
	       S.
	    */

	    if (length != nrow*nrow) {
		throwLogicError("invalid length in sampleWishart");
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
	    int ni = asInteger(nrow);
	    F77_DPOTRF("L", &ni, &C[0], &ni, &info);
	    if (info != 0) {
		throwRuntimeError("Failed to get Cholesky decomposition of R");
	    }
	    F77_DTRTRI("L", "N", &ni, &C[0], &ni, &info);
	    if (info != 0) {
		throwRuntimeError("Failed to invert Cholesky decomposition of R");
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
		Z_j[j] = sqrt(rchisq(df - j, rng));    
		for (unsigned long i = j + 1; i < nrow; i++) {
		    Z_j[i] = 0;
		}
	    }

	    // Z = Z %*% C 
	    double one = 1;
	    
	    F77_DTRMM("R", "U", "N", "N", &ni, &ni, &one, &C[0], &ni,
		      &Z[0], &ni);

	    // X = t(Z) %*% Z
	    double zero = 0;
	    F77_DSYRK("U", "T", &ni, &ni, &one, &Z[0], &ni, &zero, X, &ni);

	    // Copy upper triangle of X to lower triangle
	    for (unsigned long i = 0; i < nrow; ++i) {
		for (unsigned long j = 0; j < i; ++j) {
		    X[j * nrow + i] = X[i * nrow + j];
		}
	    }
	}


    } //namespace glm
} //namespace jags
