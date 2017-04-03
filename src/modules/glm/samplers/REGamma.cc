#include <config.h>

#include "REGamma.h"

#include <JRmath.h>

#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>

using std::vector;

namespace jags {
    namespace glm {

	REGamma::REGamma(SingletonGraphView const *tau,
			 GraphView const *eps,
			 vector<SingletonGraphView const *> const &sub_eps,
			 vector<Outcome *> const &outcomes,
			 unsigned int chain)
	    : REMethod(tau, eps, sub_eps, outcomes, chain)
	{
	}

	void REGamma::updateTau(RNG *rng)
	{
	    vector<Node const*> const &par = _tau->node()->parents();
	    double r = *par[0]->value(_chain); //shape
	    double mu = *par[1]->value(_chain); //rate (1/scale)
    
	    // Likelihood
	    vector<StochasticNode *> const &sch = _tau->stochasticChildren();
	    for (unsigned int i = 0; i < sch.size(); ++i) {
		double Y = *sch[i]->value(_chain);
		r += 0.5;
		mu += Y * Y / 2.0;
	    }

	    double x = rgamma(r, 1.0/mu, rng);
	    _tau->setValue(&x, 1, _chain);  
	}
    }
}
