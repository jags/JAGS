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
		//FIXME: We could use blas call dsyr here
		for (int j = 0; j < m; j++) {
		    for (int k = 0; k < m; k++) {
			R[j*m + k] += Y[j] * Y[k];
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
	    vector<double> A(m2, 0);
	    for (unsigned int j = 0; j < m; ++j) {
		A[j * m + j] = 1/(S[j] * S[j]);
	    }

	    vector<double> b(m, 0);
	    unsigned int N = _outcomes.size();
	    for (unsigned int i = 0; i < N; ++i) {
		double Y = _outcomes[i]->value();
		double mu = _outcomes[i]->mean();
		double lambda = _outcomes[i]->precision();
		vector<double> X(m);
		for (unsigned int j = 0; j < m; ++j) {
		    X[j] =  Zx[j*N+i]/sigma0[j];
		}
		for (unsigned int j = 0; j < m; ++j) {
		    for (unsigned int k = 0; k < m; ++k) {
			A[j*m + k] += X[j] * X[k] * lambda;
		    }
		    b[j] += (Y - mu) * X[j] * lambda;
		}
	    }
	    
	    //Sample each sigma from its full conditional
	    //Fixme: wouldn't it be better to do block sampling here?
	    //Fixme: not reversible
	    for (unsigned int j = 0; j < m; ++j) {
		double mu  = _sigma[j] + b[j]/A[j*m+j];
		double sigma = sqrt(1.0/A[j*m+j]);
		_sigma[j] = lnormal(0, rng, mu, sigma);
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
	    vector<double> eps(_eps->length());
	    _eps->getValue(eps, _chain);
	    unsigned int Neff = _eps->nodes().size();
	    for (unsigned int i = 0; i < Neff; ++i) {
		for (unsigned int j = 0; j < m; ++j) {
		    eps[m*i + j] *= sigma_ratio[j];
		}
	    }
	    _eps->setValue(eps, _chain);

	    /*
	    //Rescale tau
	    double tau = *_tau->node()->value(_chain);
	    tau /= (sigma_ratio * sigma_ratio);
	    _tau->setValue(&tau, 1, _chain);
	    */
	}

	void REScaledWishart::update(RNG *rng) {
	    
	    // Update outcomes
	    for (vector<Outcome*>::const_iterator p = _outcomes.begin();
		 p != _outcomes.end(); ++p)
	    {
		(*p)->update(rng);
	    }
	    
	    updateEps(rng);
	    updateTau(rng); //Sufficient parameterization
	    updateSigma(rng); //Ancillary parameterization
	    updateTau(rng); //Sufficient parameterization
	}

    }
}
