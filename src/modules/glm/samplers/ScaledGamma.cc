#include <config.h>
#include "ScaledGamma.h"

#include <graph/StochasticNode.h>
#include <graph/Graph.h>
#include <sampler/Linear.h>
#include <sampler/SingletonGraphView.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>

#include <vector>
#include <cmath>
#include <algorithm>

#include <JRmath.h>

using std::vector;
using std::fill;

//debuggin
#include <iostream>

namespace jags {
    namespace glm {

	static inline double 
	getScale(StochasticNode const *snode, unsigned int chain)
	{
	    return *snode->parents()[1]->value(chain);
	}

	void ScaledGamma::calCoef()
	{   
	    const double xold = _gv->node()->value(_chain)[0];
	    vector<StochasticNode *> const &sch = _gv->stochasticChildren();

	    for (unsigned int i = 0; i < sch.size(); ++i) {
		_coef[i] = getScale(sch[i], _chain);
	    }
	    double val = 2 * xold;
	    _gv->setValue(&val, 1, _chain);
	    for (unsigned int i = 0; i < sch.size(); ++i) {
		if (_coef[i] == getScale(sch[i], _chain)) {
		    _coef[i] = 0;
		}
		else {
		    _coef[i] /= xold;
		}
	    }
	    _gv->setValue(&xold, 1, _chain);
	}


	ScaledGamma::ScaledGamma(SingletonGraphView const *gv,
				 unsigned int chain)
	    : SampleMethodNoAdapt(), _gv(gv), _chain(chain),
	      _coef(gv->stochasticChildren().size())
	{
	    if(gv->deterministicChildren().empty()) {
		fill(_coef.begin(), _coef.end(), 1);
		_fast = true;
	    }
	    else if (checkScale(gv, true)) {
		//One-off calculation of fixed scale transformation
		calCoef();
		_fast = true;
	    }
	    else {
		_fast = false;
	    }
	    
	    //Initialize hyper-parameter _a to its prior mean
	    vector<Node const*> const &par = gv->node()->parents();
	    double S = *par[0]->value(chain); //Prior scale
	    double df = *par[1]->value(chain); //Prior degrees of freedom
	    
	    double x = gv->node()->value(chain)[0];
	    double a_shape = (1 + df)/2; // shape
	    double a_rate = df * x + 1/(S*S); // 1/scale
	    _a = a_shape/a_rate;

	}

	bool ScaledGamma::canSample(StochasticNode *snode, Graph const &graph)
	{
	    if (snode->distribution()->name() != "dscaled.gamma")
		return false;

	    if (isBounded(snode)) return false;
    
	    SingletonGraphView gv(snode, graph);
	    
	    // Check stochastic children
	    vector<StochasticNode *> const &stoch_nodes = 
		gv.stochasticChildren();
	    for (unsigned int i = 0; i < stoch_nodes.size(); ++i) {
		if (isBounded(stoch_nodes[i])) {
		    return false; 
		}
		if (stoch_nodes[i]->distribution()->name() != "dnorm") {
		    return false;
		}
		if (gv.isDependent(stoch_nodes[i]->parents()[0])) {
		    return false; //non-scale parameter depends on snode
		}
	    }
	    
	    // Check deterministic descendants are scale transformations 
	    if (!checkScale(&gv, false)) {
		return false;
	    }
	    return true; //We made it!
	}

	static void
	sample_gamma(double &x, double shape, double rate, RNG *rng,
	    bool relax)
	{
	    double scale = 1/rate;

	    if (relax) {
		double K = 19;
		
		double u = pgamma(x, shape, scale, true, false);
		double r = rbinom(K, u, rng);
		if (r > K - r) {
		    u *= rbeta(K - r + 1, 2*r -  K, rng);
		}
		else if (r < K - r) {
		    u = 1 - (1-u) * rbeta(r + 1, K - 2*r, rng);
		}
		x = qgamma(u, shape, scale, true, false);
	    }
	    else {
		x = rgamma(shape, scale, rng);
	    }
	}

	void ScaledGamma::update(RNG *rng)
	{
	    vector<StochasticNode *> const &sch = _gv->stochasticChildren();
	    unsigned int nchildren = sch.size();
    
	    double r = 0; //shape
	    double mu = 0; //rate (1/scale)
    
	    // Likelihood

	    if (!_fast) calCoef();
		
	    for (unsigned int i = 0; i < nchildren; ++i) {

		if (_coef[i] == 0) continue;
		
		StochasticNode const *schild = sch[i];
		vector<Node const*> const &cparam = schild->parents();
		double Y = *schild->value(_chain);
		double m = *cparam[0]->value(_chain); //location parameter 

		r += 0.5;
		mu += _coef[i] * (Y - m) * (Y - m) / 2;
	    }

	    vector<Node const*> const &par = _gv->node()->parents();
	    double S = *par[0]->value(_chain);
	    double df = *par[1]->value(_chain);

	    double x = _gv->node()->value(_chain)[0];

	    bool first = rng->uniform() < 0.5 ;
	    sample_gamma(_a, (1 + df)/2, df * x + 1/(S*S), rng, first);
	    sample_gamma( x, r + df/2, mu + df * _a, rng, true);
	    sample_gamma(_a, (1 + df)/2, df * x + 1/(S*S), rng, !first);

	    _gv->setValue(&x, 1, _chain);  
	}

    }
}
