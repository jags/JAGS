#include <config.h>

#include "REScaledGamma2.h"

#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <rng/TruncatedNormal.h>
#include <JRmath.h>

#include <cmath>

using std::vector;
using std::sqrt;

extern cholmod_common *glm_wk;

namespace jags {
    namespace glm {

	REScaledGamma2::REScaledGamma2(SingletonGraphView const *tau,
				       GLMMethod const *glmmethod)
	    : REMethod2(tau, glmmethod)
	{
	    //Initialize hyper-parameter _sigma 
	    vector<Node const*> const &par = tau->node()->parents();
	    double S = *par[0]->value(_chain); //Prior scale
	    double df = *par[1]->value(_chain); //Prior degrees of freedom

	    double x = tau->node()->value(_chain)[0];
	    double a_shape = (1 + df)/2; // shape
	    double a_rate = df * x + 1/(S*S); // 1/scale
	    _sigma = sqrt(a_shape/a_rate);
	}

	void REScaledGamma2::updateTau(RNG *rng)
	{
	    double df = *_tau->node()->parents()[1]->value(_chain);

	    // Prior
	    double shape = df/2.0; 
	    double rate = df * _sigma * _sigma / 2.0; // 1/scale
    
	    // Likelihood
	    vector<StochasticNode *> const &eps = _tau->stochasticChildren();
	    for (unsigned int i = 0; i < eps.size(); ++i) {
		double Y = *eps[i]->value(_chain);
		double mu = *eps[i]->parents()[0]->value(_chain);
		shape += 0.5;
		rate += (Y - mu) * (Y - mu) / 2.0;
	    }
	    
	    double x = rgamma(shape, 1.0/rate, rng);
	    _tau->setValue(&x, 1, _chain);  
	}

	void REScaledGamma2::updateSigma(RNG *rng)
	{
	    double sigma0 = _sigma;

	    calDesignSigma();

	    //Prior scale
	    vector<Node const*> const &par = _tau->node()->parents();
	    double S = *par[0]->value(_chain);

	    double const *Zx = static_cast<double const *>(_z->x);
	    
	    //Get parameters of posterior distribution for _sigma
	    //Precision is A and mean is b/A
	    double priorprec = 1.0/(S*S);
	    double A = priorprec;
	    double b = - _sigma * priorprec;

	    calCoefSigma(&A, &b, &_sigma, 1);

	    //Set new value of sigma
	    _sigma = lnormal(0, rng, _sigma + b/A, 1/sqrt(A));

	    //Rescale tau
	    double tau = *_tau->node()->value(_chain);
	    tau *= (sigma0 * sigma0)/(_sigma * _sigma);
	    _tau->setValue(&tau, 1, _chain);
	}
	
	bool REScaledGamma2::isAdaptive() const
	{
	    return false;
	}
	
	void REScaledGamma2::adaptOff()
	{
	}
	
	bool REScaledGamma2::checkAdaptation() const
	{
	    return true;
	}
	
    }
}
