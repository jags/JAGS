#include "REGammaSlicer.h"
#include "Outcome.h"

#include <JRmath.h>

using std::copy;

namespace jags {
    namespace glm {
	
	REGammaSlicer::REGammaSlicer(vector<Outcome*> const &outcomes,
				     cholmod_sparse const *x,
				     cholmod_dense const *z,
				     double const *shape,
				     double const *rate,
				     double sigma)
	    : Slicer(1, 10), _outcomes(outcomes.size()), _x(x), _z(z),
	      _shape(shape), _rate(rate), _sigma(sigma), _sigma0(sigma)
	{
	    copy(outcomes.begin(), outcomes.end(), _outcomes.begin());
	}
	
	
	double REGammaSlicer::value() const
	{
	    return _sigma;
	}
	void REGammaSlicer::setValue(double x)
	{
	    _sigma = x;
	}

	void REGammaSlicer::setSigma(double x)
	{
	    _sigma0 = _sigma = x;
	}
	
	void REGammaSlicer::getLimits(double *lower, double *upper) const
	{
	    *lower = 1e-06;
	    *upper = 1e06;
	}
	
	double REGammaSlicer::logDensity() const
	{
	    double const *Zx = static_cast<double const *>(_z->x);

	    double tau = 1/(_sigma * _sigma);
	    double scale = 1/(*_rate);
	    double logprior = dgamma(tau, *_shape, scale, 1) - 3 * log(_sigma);
	    
	    double loglik = 0;
	    unsigned int N = _outcomes.size();
	    for (unsigned int i = 0; i < N; ++i) {
		double Y = _outcomes[i]->value();
		double mu = _outcomes[i]->mean();
		double lambda = _outcomes[i]->precision();
		double delta = Y - mu - (_sigma/_sigma0 - 1) * Zx[i];
		loglik -= lambda * delta * delta / 2;
	    }
	    return logprior + loglik;
	}

	void REGammaSlicer::update(RNG *rng)
	{
	    updateStep(rng);
	}
	
    }
}
       
