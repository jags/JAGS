#include <config.h>

#include "REGamma2.h"

#include <JRmath.h>
#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>

#include <cmath>

using std::vector;
using std::sqrt;


namespace jags {

    //Utility functions used by the constructor

    static inline double const *
    SHAPE(SingletonGraphView const *tau, unsigned int chain)
    {
	return tau->node()->parents()[0]->value(chain);
    }

    static inline double const *
    RATE(SingletonGraphView const *tau, unsigned int chain)
    {
	return tau->node()->parents()[1]->value(chain);
    }

    static inline double 
    SIGMA(SingletonGraphView const *tau, unsigned int chain)
    {
	return 1.0/sqrt(*tau->node()->value(chain));
    }

    namespace glm {

	class Outcome;
	
	REGamma2::REGamma2(SingletonGraphView const *tau,
			   GLMMethod const *glmmethod)
	    : REMethod2(tau, glmmethod),
	      _slicer(this, SHAPE(tau, _chain), RATE(tau, _chain),
		      SIGMA(tau, _chain))
	{
	}

	void REGamma2::updateTau(RNG *rng)
	{
	    vector<Node const*> const &par = _tau->node()->parents();
	    double shape = *par[0]->value(_chain); 
	    double rate = *par[1]->value(_chain); //(1/scale)

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

	void REGamma2::updateSigma(RNG *rng)
	{
	    double tau = _tau->node()->value(_chain)[0];
	    double sigma0 = 1/sqrt(tau);

	    calDesignSigma();
	    
	    _slicer.setSigma(sigma0);
	    _slicer.update(rng);
	    double sigma1 = _slicer.value();

	    //Set new value of precision parameter
	    double x = 1/(sigma1 * sigma1);
	    _tau->setValue(&x, 1, _chain);
	}

	bool REGamma2::isAdaptive() const
	{
	    return true;
	}
	
	void REGamma2::adaptOff()
	{
	    _slicer.adaptOff();
	}
	
	bool REGamma2::checkAdaptation() const
	{
	    return _slicer.checkAdaptation();
	}

    }
}
