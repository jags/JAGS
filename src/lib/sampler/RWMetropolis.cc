#include <config.h>
#include <sampler/RWMetropolis.h>
#include <rng/RNG.h>

#include <cmath>

using std::vector;
using std::log;
using std::exp;
using std::fabs;
using std::isfinite;
using std::min;

namespace jags {
    
    /* FIXME: Utility functions copied from src/modules/bugs/samplers/MNormal.cc */
    
    // Target effective sample size 
    static double ESS(unsigned long t, double b, double n0) {
	if (t <= 1) {
	    return n0 + t;
	} else {
	    return n0 + 1 + b*(t-1);
	}
    }
    
    // Returns learning rate for the sample mean and variance for
    // a given sample size t.
    static double solve_lambda(unsigned long t, double b, double n0) {
	double S0 = ESS(t-1, b, n0);
	double S1 = ESS(t, b, n0);
	double delta = S1 - S0;
	
	return (1 + sqrt(S0*(1-delta)/S1))/(1 + S0);
    }


RWMetropolis::RWMetropolis(vector<double> const &value, double step,
			   double a, double delta, double nstart, double min_step)
    : Metropolis(value), _step_adapter(step, a, delta, nstart, min_step), _niter(0), _pmean(0)
{
}

RWMetropolis::~RWMetropolis()
{
}

void RWMetropolis::rescale(double p)
{
    p = min(p, 1.0);
    _step_adapter.rescale(p);

    _niter++;
    _pmean += solve_lambda(_niter, 0.5, 100) * (p - _pmean);
}

void RWMetropolis::update(RNG *rng)
{
    vector<double> value(length());
    getValue(value);

    double log_p = logDensity() + logJacobian(value);
    step(value, _step_adapter.stepSize(), rng);
    setValue(value);
    double log_p_new = logDensity() + logJacobian(value);
    double odds = ( isfinite( log_p ) && isfinite( log_p_new ) )
                ?  exp( log_p_new - log_p )
                : ( log_p_new > log_p ? 1 : 0 );
    accept(rng, odds);
}

bool RWMetropolis::checkAdaptation() const
{
    if (_niter < 100) return false;

    return fabs(_step_adapter.logitDeviation(_pmean)) < 0.5;
}

void RWMetropolis::step(vector<double> &value, double s, RNG *rng) const
{
    for (unsigned int i = 0; i < value.size(); ++i) {
	value[i] += rng->normal() * s;
    }
}

double RWMetropolis::logJacobian(vector<double> const &) const
{
    return 0;
}

} //namespace jags
