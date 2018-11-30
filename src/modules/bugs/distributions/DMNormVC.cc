#include <config.h>
#include <util/dim.h>
#include <util/nainf.h>
#include "DMNormVC.h"
#include "DMNorm.h"

#include <lapack.h>
#include <matrix.h>

#include <cmath>
#include <vector>
#include <cfloat>

#include <JRmath.h>

#include "matrix.h"

using std::vector;

namespace jags {
    namespace bugs {
    
	DMNormVC::DMNormVC()
	    : ArrayDist("dmnorm.vcov", 2) 
	{}

	double
	DMNormVC::logDensity(double const *x, PDFType type,
			     vector<double const *> const &parameters,
			     vector<vector<unsigned long> > const &dims) const
	{
	    double const * mu = parameters[0];
	    double const * V  = parameters[1];
	    unsigned long m = dims[0][0];
	    
	    vector<double> T(m * m);
	    inverse_chol (&T[0], V, m);

	    double loglik = 0;
	    vector<double> delta(m);
	    for (unsigned long i = 0; i < m; ++i) {
		delta[i] = x[i] - mu[i];
		loglik -= (delta[i] * T[i + i * m] * delta[i])/2;
		for (unsigned long j = 0; j < i; ++j) {
		    loglik -= delta[i] * T[i + j * m] * delta[j];
		}
	    }

	    switch(type) {
	    case PDF_PRIOR:
		break;
	    case PDF_LIKELIHOOD:
		loglik -= logdet(V, m)/2;
		break;
	    case PDF_FULL:
		loglik -= logdet(V, m)/2 + m * M_LN_SQRT_2PI;
		break;
	    }
    
	    return loglik;
	}

	void
	DMNormVC::randomSample(double *x,
			       vector<double const *> const &parameters,
			       vector<vector<unsigned long> > const &dims,
			       RNG *rng) const
	{
	    double const * mu = parameters[0];
	    double const * T = parameters[1];
	    unsigned long m = dims[0][0];
	    
	    DMNorm::randomsample(x, mu, T, false, m, rng);
	}

	bool
	DMNormVC::checkParameterDim(vector<vector<unsigned long> > const &dims)
	    const
	{
	    //Allow scalar mean and precision. 
	    if (isScalar(dims[0]) && isScalar(dims[1])) return true;

	    //Vector mean and matrix precision
	    if (!isVector(dims[0])) return false;
	    if (!isSquareMatrix(dims[1])) return false;
	    if (dims[0][0] != dims[1][0]) return false;
    
	    return true;
	}

	vector<unsigned long>
	DMNormVC::dim(vector<vector<unsigned long> > const &dims) const
	{
	    return dims[0];
	}
	
	void
	DMNormVC::support(double *lower, double *upper,
			  vector<double const *> const &,
			  vector<vector<unsigned long> > const &dims) const
	{
	    unsigned length = dims[0][0];
	    for (unsigned long i = 0; i < length; ++i) {
		lower[i] = JAGS_NEGINF;
		upper[i] = JAGS_POSINF;
	    }
	}

	bool DMNormVC::isSupportFixed(vector<bool> const &) const
	{
	    return true;
	}

	bool DMNormVC::checkParameterValue(vector<double const *> const &,
					   vector<vector<unsigned long> > const &) const
	{
	    return true; //FIXME: define in base class
	}

	void DMNormVC::randomsample(double *x, double const *mu, double const *Sigma,
				    unsigned long nrow,
				    vector<bool> const &observed, unsigned long nobs,
				    RNG *rng)
	{
	    if (nrow == nobs) return;
	    
	    unsigned long nfree = nrow - nobs;
	    unsigned long N = nfree*nfree;
	    
	    vector<double> Sigff(nfree*nfree), Sigfo(nfree*nobs), Sigof(nobs*nfree), Sigoo(nobs*nobs);
	    vector<double> wf(nfree), wo(nobs);
		
	    if (nobs == 0) {
		copy(Sigma, Sigma + N, Sigff.begin());
	    }
	    else {
		//Partition variance matrix Sigma into free and observed parts
		for (unsigned long i = 0, pf = 0, po = 0; i < nrow; ++i) {
		    if (observed[i]) {
			for (unsigned long j = 0, qf = 0, qo = 0; j < nrow; ++j) {
			    if (observed[j]) {
				Sigoo[po*nobs + qo++] = Sigma[i*nrow + j];
			    }
			    else {
				Sigof[po*nobs + qf++] = Sigma[i*nrow + j];
			    }
			}
			wo[po++] = x[i] - mu[i];
		    }
		    else {
			for (unsigned long j = 0, qf = 0, qo = 0; j < nrow; ++j) {
			    if (observed[j]) {
				Sigfo[pf*nfree + qo++] = Sigma[i*nrow + j];
			    }
			    else {
				Sigff[pf*nfree + qf++] = Sigma[i*nrow + j];
			    }
			}
			wf[pf++] = mu[i];
		    }
		}
	    }

	    //Factorize Sigff.
	    //...

	    if (nobs > 0) {
		//Solve Sigma[o,o] %*% x = (Y[o] - mu[o])
		//Do matrix multiplication Sigma[f,o] %*% x + mu[f]
	    }
		
	    //Generate independent normal random variables eps
	    //Calculate L %*% eps to get correct variance

	    //Write back to free portions of x
	    for (unsigned long i = 0, p = 0; i < nrow; ++i) {
		if (!observed[i]) x[i] = wf[p++];
	    }
			
		    
		    
	}

	
    } //namespace bugs
} //namespace jags
