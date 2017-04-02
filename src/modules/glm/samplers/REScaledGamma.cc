#include <config.h>

#include "REScaledGamma.h"
#include "Outcome.h"
#include "RESampler.h"

#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
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
	    //Initialize hyper-parameter _a to its prior mean
	    vector<Node const*> const &par = tau->node()->parents();
	    double S = *par[0]->value(chain); //Prior scale
	    double df = *par[1]->value(chain); //Prior degrees of freedom

	    double x = tau->node()->value(chain)[0];
	    double a_shape = (1 + df)/2; // shape
	    double a_rate = df * x + 1/(S*S); // 1/scale
	    _a = a_shape/a_rate;
	}

	static inline double sample_gamma(double shape, double rate, RNG *rng)
	{
	    return rgamma(shape, 1/rate, rng);
	}

	void REScaledGamma::updateTau(RNG *rng)
	{
	    double r = 0; //shape
	    double mu = 0; //rate (1/scale)
    
	    // Likelihood
	    vector<StochasticNode *> const &sch = _tau->stochasticChildren();
	    for (unsigned int i = 0; i < sch.size(); ++i) {
		double Y = *sch[i]->value(_chain);
		r += 0.5;
		mu += Y * Y / 2.0;
	    }

	    // Prior
	    vector<Node const*> const &par = _tau->node()->parents();
	    double S = *par[0]->value(_chain);
	    double df = *par[1]->value(_chain);

	    double x = _tau->node()->value(_chain)[0];

	    _a = sample_gamma((1 + df)/2.0, df * x + 1/(S*S), rng);
	    x = sample_gamma(r + df/2.0, mu + df * _a, rng);
	    _a = sample_gamma((1 + df)/2.0, df * x + 1/(S*S), rng);

	    _tau->setValue(&x, 1, _chain);  
	}
    }
}
