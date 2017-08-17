#include "REGammaSlicer2.h"
#include "REGamma2.h"

#include <JRmath.h>

namespace jags {
    namespace glm {
	
	REGammaSlicer2::REGammaSlicer2(REGamma2 const *regamma,
				     double const *shape,
				     double const *rate,
				     double sigma)
	    : Slicer(1, 10), _regamma(regamma),
	      _shape(shape), _rate(rate), _sigma(sigma), _sigma0(sigma)
	{
	}
	
	double REGammaSlicer2::value() const
	{
	    return _sigma;
	}
	void REGammaSlicer2::setValue(double x)
	{
	    _sigma = x;
	}

	void REGammaSlicer2::setSigma(double x)
	{
	    _sigma0 = _sigma = x;
	}
	
	void REGammaSlicer2::getLimits(double *lower, double *upper) const
	{
	    *lower = 1e-06;
	    *upper = 1e06;
	}
	
	double REGammaSlicer2::logDensity() const
	{
	    /* Log gamma prior on tau, the precision parameter, with a
	       Jacobian term for the transformation to the standard
	       deviation parameter _sigma */
	    double tau = 1/(_sigma * _sigma);
	    double scale = 1/(*_rate);
	    double logprior = dgamma(tau, *_shape, scale, 1) - 3 * log(_sigma);

	    /* For the likelihood, We refer back to the REGamma object
	       that owns this REGammaSlicer2 object. */
	    double loglik = _regamma->logLikelihoodSigma(&_sigma, &_sigma0, 1);

	    return logprior + loglik;
	}

	void REGammaSlicer2::update(RNG *rng)
	{
	    updateStep(rng);
	}
	
    }
}
       
