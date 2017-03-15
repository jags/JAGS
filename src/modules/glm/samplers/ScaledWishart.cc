#include <config.h>

#include "ScaledWishart.h"

#include <rng/RNG.h>
#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <graph/MixtureNode.h>
#include <sampler/Linear.h>
#include <sampler/SingletonGraphView.h>
#include <module/ModuleError.h>

#include <set>
#include <vector>
#include <cmath>
#include <algorithm>

#include <JRmath.h>

using std::vector;
using std::set;
using std::sqrt;
using std::string;
using std::copy;

#define F77_DPOTRF F77_FUNC(dpotrf,DPOTRF)
#define F77_DPOTRI F77_FUNC(dpotri, DPOTRI)

extern "C" {
    void F77_DPOTRF (const char *uplo, const int *n, double *a,
		     const int *lda, const int *info);
    void F77_DPOTRI (const char *uplo, const int *n, double *a,
		     const int *lda, const int *info);
}


static bool inverse_spd (double *A, int n)
{
    /* invert n x n symmetric positive definite matrix A*/

    int info = 0;
    F77_DPOTRF ("L", &n, A, &n, &info);
    if (info != 0) {
	return false;
    }

    F77_DPOTRI ("L", &n, A, &n, &info); 

    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < i; ++j) {
	    A[i*n + j] = A[j*n + i];
	}
    }

    return info == 0;
}

//FIXME. Bah! This is copy-pasted from the bugs module but it is
//clearly inefficient. We invert the matrix using the Cholesky
//decomposition and then get the Cholesky decomposition of the
//inverse!  This whole function should be done with BLAS and LAPACK
//calls.

static void sampleWishart(double *x, int length,
			  double const *R, double k, int nrow,
			  jags::RNG *rng)
{
    /* 
       Generate random Wishart variable, using an algorithm proposed
       by Bill Venables and originally implemented in S.
    */

    if (length != nrow*nrow) {
	jags::throwLogicError("invalid length in DWish::randomSample");
    }

    /* 
       Get inverse of R. Venables' algorithm was implemented in
       terms of the inverse of R, but we use a different parameterization
       to preserve conjugacy.
    */
    double * C = new double[length];
    copy(R, R + length, C);
    if(!inverse_spd(C, nrow)) {
	jags::throwRuntimeError("Inverse failed in DWish::randomSample");
    }
    /* Get Choleskly decomposition of C */
    int info = 0;
    F77_DPOTRF("U", &nrow, C, &nrow, &info);
    if (info != 0) {
	jags::throwRuntimeError("Failed to get Cholesky decomposition of R");
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


namespace jags {
    namespace glm {

	/*
	static inline double getPrecision0(StochasticNode const *snode, 
					   unsigned int chain)
	{
	    //Returns the first element of the precision matrix for a node
	    //with a multivariate normal distribution.
	    
	    return snode->parents()[1]->value(chain)[0];
	}
	*/

	bool
	ScaledWishart::canSample(StochasticNode *snode, Graph const &graph)
	{
	    if (snode->distribution()->name() != "dscaled.wishart")
		return false;

	    if (isBounded(snode))
		return false;
  
	    SingletonGraphView gv(snode, graph);
	    vector<StochasticNode *> const &schild = gv.stochasticChildren();

	    // Check stochastic children
	    for (unsigned int i = 0; i < schild.size(); ++i) {
		if (isBounded(schild[i])) {
		    return false; //Bounded
		}
		if (schild[i]->distribution()->name() != "dmnorm" &&
		    schild[i]->distribution()->name() != "dnorm")
		{
		    return false;
		}
		if (schild[i]->parents()[1] != snode) {
		    return false;
		}
		if (gv.isDependent(schild[i]->parents()[0])) {
		    return false; //mean parameter depends on snode
		}
	    }

	    if (!gv.deterministicChildren().empty()) {
		return false;
	    }

	    return true;
	}

	ScaledWishart::ScaledWishart(SingletonGraphView const *gv,
					   unsigned int chain)
	    : SampleMethodNoAdapt(), _gv(gv), _chain(chain)
	{
	    vector<Node const *> const &par = _gv->node()->parents();  
	    unsigned int nrow = par[0]->dim()[0];
	    
	    double const *S = par[0]->value(_chain); //Prior scale
	    double df = *par[1]->value(_chain); //Prior degrees of freedom
	    double const *x = gv->node()->value(_chain);

	    _a = vector<double>(nrow);	    
	    for (unsigned int k = 0; k < nrow; ++k) {
		double shape = (1 + df)/2; // shape
		double rate = df * x[k*nrow+k] + 1/(S[k]*S[k]); // 1/scale
		_a[k] = shape/rate;
	    }
	}
	
	void ScaledWishart::update(RNG *rng)
	{
	    vector<Node const *> const &param = _gv->node()->parents();  

	    double tdf = *param[1]->value(_chain);
	    double const *S = param[0]->value(_chain);
	    int nrow = param[0]->dim()[0];
	    int N = nrow * nrow;
	    double const *x = _gv->node()->value(_chain);
	    
	    //Update hyper-parameters _a
	    for (int i = 0; i < nrow; ++i) {
		double shape = (nrow + tdf)/2;
		double rate = tdf * x[i*nrow+i] + 1/(S[i]*S[i]);
		_a[i] = rgamma(shape, 1/rate, rng);
	    }

	    //Prior contribution to Wishart parameters
	    double wdf = nrow + tdf - 1; //Degrees of freedom for Wishart
	    vector<double> R(N, 0); //Scale matrix for Wishart
	    for (int i = 0; i < nrow; ++i) {
		R[i*nrow+i] = 2 * tdf * _a[i];
	    }

	    //Likelihood contribution to Wishart parameters
	    vector<StochasticNode *> const &schild = _gv->stochasticChildren();

	    for (vector<StochasticNode*>::const_iterator p = schild.begin();
		 p != schild.end(); ++p)
	    {
		double const *Y = (*p)->value(_chain);
		double const *mu = (*p)->parents()[0]->value(_chain);
	    
		for (int i = 0; i < nrow; i++) {
		    for (int j = 0; j < nrow; j++) {
			R[i*nrow + j] += (Y[i] - mu[i]) * (Y[j] - mu[j]);
		    }
		}
		wdf += 1;
	    }

	    vector<double> xnew(N);
	    sampleWishart(&xnew[0], N, &R[0], wdf, nrow, rng);
	    _gv->setValue(xnew, _chain);
	}

    }
}
