#include <config.h>

#include "REScaledWishart2.h"
#include "Outcome.h"
#include "RESampler.h"
#include "SampleWishart.h"

#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <rng/TruncatedNormal.h>
#include <JRmath.h>

#include <cmath>

using std::vector;
using std::sqrt;

extern cholmod_common *glm_wk;

namespace jags {
    namespace glm {

	REScaledWishart2::REScaledWishart2(SingletonGraphView const *tau,
					   GLMMethod const *glmmethod)
	    : REMethod2(tau, glmmethod)
	{
	    vector<Node const*> const &par = tau->node()->parents();
	    double const *S = par[0]->value(_chain); //Prior scale
	    unsigned int nrow = par[0]->length();
	    double tdf = *par[1]->value(_chain); //Prior degrees of freedom
	    double const *x = tau->node()->value(_chain);

	    //Initialize hyper-parameter _sigma
	    _sigma = vector<double>(nrow);
	    for (unsigned int j = 0; j < nrow; ++j) {
		double a_shape = (nrow + tdf)/2.0;
		double a_rate = tdf * x[j + nrow*j] + 1.0/(S[j]*S[j]);
		double a = a_shape/a_rate; 
		_sigma[j] = sqrt(2*a);
	    }
	}

	void REScaledWishart2::updateTau(RNG *rng)
	{
	    int m = _sigma.size();
	    int m2 = m * m;
	    double tdf = *_tau->node()->parents()[1]->value(_chain);

	    //Prior 
	    double wdf = m + tdf - 1; //Degrees of freedom for Wishart
	    vector<double> R(m2, 0); //Scale matrix for Wishart
	    for (int j = 0; j < m; ++j) {
		R[j*m + j] = tdf * _sigma[j] * _sigma[j];
	    }

	    //Likelihood
	    vector<StochasticNode *> const &eps = _tau->stochasticChildren();
	    for (vector<StochasticNode*>::const_iterator p = eps.begin();
		 p != eps.end(); ++p)
	    {
		double const *Y = (*p)->value(_chain);
		double const *mu = (*p)->parents()[0]->value(_chain);
		//FIXME: We could use blas call dsyr here
		for (int j = 0; j < m; j++) {
		    for (int k = 0; k < m; k++) {
			R[j*m + k] += (Y[j] - mu[j]) * (Y[k] - mu[k]);
		    }
		}
		wdf += 1;
	    }

	    vector<double> xnew(m2);
	    sampleWishart(&xnew[0], m2, &R[0], wdf, m, rng);
	    _tau->setValue(xnew, _chain);
	}

	void REScaledWishart2::updateSigma(RNG *rng)
	{
	    vector<double> sigma0 = _sigma;
	    calDesignSigma();

	    //Prior scale
	    vector<Node const*> const &par = _tau->node()->parents();
	    double const *S = par[0]->value(_chain);

	    unsigned int m  = _z->ncol;
	    unsigned int m2 = m * m;
	    
	    //Get parameters of posterior distribution for _sigma
	    //Precision is A and mean is inverse(A) %*% b
	    //We work on a scale where the current value sigma0 is the origin
	    vector<double> A(m2, 0);
	    vector<double> b(m, 0);
	    for (unsigned int j = 0; j < m; ++j) {
		double priorprec = 1.0/(S[j] * S[j]);
		A[j * m + j] = priorprec;
		b[j] = - sigma0[j] * priorprec;
	    }

	    calCoefSigma(&A[0], &b[0], &sigma0[0], m);
	    
	    //Sample each sigma from its full conditional
	    //Fixme: wouldn't it be better to do block sampling here?
	    //Fixme: not reversible
	    for (unsigned int j = 0; j < m; ++j) {
		double sigma_mean  = _sigma[j] + b[j]/A[j*m+j];
		double sigma_sd = sqrt(1.0/A[j*m+j]);
		_sigma[j] = lnormal(0, rng, sigma_mean, sigma_sd);
		double delta = _sigma[j] - sigma0[j];
		for (unsigned int k = 0; k < m; ++k) {
		    b[k] -= delta * A[m*j + k];
		}
	    }

	    //Rescale tau
	    double const *tau0 = _tau->node()->value(_chain);
	    vector<double> scale(m);
	    for (unsigned int j = 0; j < m; ++j) {
		scale[j] = sigma0[j]/_sigma[j];
	    }
	    
	    vector<double> tau_scaled(m2);
	    for (unsigned int j = 0; j < m; ++j) {
		for (unsigned int k = 0; k < m; ++k) {
		    tau_scaled[j*m+k] = tau0[j*m+k] * scale[j] * scale[k];
		}
	    }
	    _tau->setValue(tau_scaled, _chain);
	}

	bool REScaledWishart2::isAdaptive() const
	{
	    return false;
	}
	
	void REScaledWishart2::adaptOff()
	{
	}
	
	bool REScaledWishart2::checkAdaptation() const
	{
	    return true;
	}
    }
}
