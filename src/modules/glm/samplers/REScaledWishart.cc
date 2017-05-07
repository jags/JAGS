#include <config.h>

#include "REScaledWishart.h"
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

	REScaledWishart::REScaledWishart(SingletonGraphView const *tau,
	    GraphView const *eps,
	    vector<SingletonGraphView const *> const &sub_eps,
	    vector<Outcome *> const &outcomes,
	    unsigned int chain)
	    : REMethod(tau, eps, sub_eps, outcomes, chain),
	      _sigma(eps->nodes()[0]->length())
	{
	    vector<Node const*> const &par = tau->node()->parents();
	    double const *S = par[0]->value(chain); //Prior scale
	    double tdf = *par[1]->value(chain); //Prior degrees of freedom
	    double const *x = tau->node()->value(chain);
	    //Initialize hyper-parameter _sigma
	    unsigned int m = _sigma.size();
	    for (unsigned int j = 0; j < m; ++j) {
		double a_shape = (m + tdf)/2.0;
		double a_rate = tdf * x[j + m*j] + 1.0/(S[j]*S[j]);
		double a = a_shape/a_rate; 
		_sigma[j] = sqrt(2*a);
	    }
	}

	void REScaledWishart::updateTau(RNG *rng)
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
	    vector<StochasticNode *> const &eps = _eps->nodes();
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

	void REScaledWishart::updateSigma(RNG *rng)
	{
	    vector<double> sigma0 = _sigma;
	    calDesignSigma();

	    //Prior scale
	    vector<Node const*> const &par = _tau->node()->parents();
	    double const *S = par[0]->value(_chain);

	    double const *Zx = static_cast<double const *>(_z->x);
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
		for (int k = 0; k < m; ++k) {
		    b[k] -= delta * A[m*j + k];
		}
	    }

	    vector<double> sigma_ratio(m);
	    for (unsigned int j = 0; j < m; ++j) {
		sigma_ratio[j] = _sigma[j]/sigma0[j];
	    }

	    //Rescale random effects
	    vector<StochasticNode *> const &eps = _eps->nodes();
	    vector<double> eval(_eps->length());
	    for (unsigned int i = 0; i < eps.size(); ++i) {
		double const *Y = eps[i]->value(_chain);
		double const *mu = eps[i]->parents()[0]->value(_chain);
		for (unsigned int j = 0; j < m; ++j) {
		    eval[m*i + j] = mu[j] + (Y[j] - mu[j]) * sigma_ratio[j];
		}
	    }
	    _eps->setValue(eval, _chain);
	    
	    /*
	    //Rescale tau
	    double tau = *_tau->node()->value(_chain);
	    tau /= (sigma_ratio * sigma_ratio);
	    _tau->setValue(&tau, 1, _chain);
	    */
	}
    }
}
