#include <config.h>

#include "REScaledGamma.h"
#include "Outcome.h"
#include "RESampler.h"

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

	REScaledGamma::REScaledGamma(SingletonGraphView const *tau,
	    GraphView const *eps,
	    vector<SingletonGraphView const *> const &sub_eps,
	    vector<Outcome *> const &outcomes,
	    unsigned int chain)
	    : REMethod(tau, eps, sub_eps, outcomes, chain)
	{
	    //Initialize hyper-parameter _sigma 
	    vector<Node const*> const &par = tau->node()->parents();
	    double S = *par[0]->value(chain); //Prior scale
	    double df = *par[1]->value(chain); //Prior degrees of freedom

	    double x = tau->node()->value(chain)[0];
	    double a_shape = (1 + df)/2; // shape
	    double a_rate = df * x + 1/(S*S); // 1/scale
	    _sigma = sqrt(a_shape/a_rate);
	}

	void REScaledGamma::updateTau(RNG *rng)
	{
	    double df = *_tau->node()->parents()[1]->value(_chain);

	    // Prior
	    double shape = df/2.0; 
	    double rate = df * _sigma * _sigma / 2.0; // 1/scale
    
	    // Likelihood
	    vector<StochasticNode *> const &eps = _eps->nodes();
	    for (unsigned int i = 0; i < eps.size(); ++i) {
		double Y = *eps[i]->value(_chain);
		shape += 0.5;
		rate += Y * Y / 2.0;
	    }
	    
	    double x = rgamma(shape, 1.0/rate, rng);
	    _tau->setValue(&x, 1, _chain);  
	}

	void REScaledGamma::updateSigma(RNG *rng)
	{
	    double sigma0 = _sigma;

	    calDesignSigma();

	    //Prior scale
	    vector<Node const*> const &par = _tau->node()->parents();
	    double S = *par[0]->value(_chain);

	    double const *Zx = static_cast<double const *>(_z->x);
	    
	    //Get parameters of posterior distribution for _sigma
	    //Precision is A and mean is b/A
	    double A = 1/(S*S);
	    double b = 0;
	    unsigned int N = _outcomes.size();
	    for (unsigned int i = 0; i < N; ++i) {
		double Y = _outcomes[i]->value();
		double mu = _outcomes[i]->mean();
		double lambda = _outcomes[i]->precision();

		double X =  Zx[i]/sigma0;

		A += X * X * lambda;
		b += (Y - mu + Zx[i]) * X * lambda;
	    }

	    //Set new value of sigma
	    //FIXME: Truncate or not?
	    //_sigma = rnorm(b/A, 1/sqrt(A), rng);
	    _sigma = lnormal(0, rng, b/A, 1/sqrt(A));
	    double sigma_ratio = _sigma/sigma0;
	    
	    //Rescale random effects
	    vector<double> eps(_eps->length());
	    _eps->getValue(eps, _chain);
	    for (unsigned int i = 0; i < eps.size(); ++i) {
		eps[i] *= sigma_ratio;
	    }
	    _eps->setValue(eps, _chain);

	    /*
	    //Rescale tau
	    double tau = *_tau->node()->value(_chain);
	    tau /= (sigma_ratio * sigma_ratio);
	    _tau->setValue(&tau, 1, _chain);
	    */
	}

	void REScaledGamma::update(RNG *rng) {
	    
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
