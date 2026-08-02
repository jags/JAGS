#include <config.h>
#include <sampler/StepAdapter.h>
#include <JRmath.h>

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <string>

using std::exp;
using std::log;
using std::logic_error;
using std::max;
using std::min;

namespace jags {




    StepAdapter::StepAdapter(double step, double a, double delta, double nstart, double min_step)
	: _a(a), _delta(delta), _theta0(log(step)), _theta(_theta0), _min_step(min_step),
	  _n(0), _nstart0(nstart), _nstart(nstart)
    {
	if (step <= 0 || a <= 0 || a >= 1 || delta <= 0 || nstart <= 0)
	    throw logic_error("Invalid initial values in StepAdapter");
    }
    
    void StepAdapter::rescale(double p)
    {
	/* Rescale step size using Robbins-Munro (1951) stochastic search algorithm to
	   reach optimal acceptance probability */
	
	p = min(p, 1.0);

	_theta += _delta * (p - _a)/(_n + _nstart);
	if (_min_step > 0) {
	    // If a minimum step size is specified, ensure we do not go below it (Spencer 2021)
	    _theta = max(_theta, log(_min_step));
	}
	
	if (abs(_theta0 - _theta) > log(3.0)) {
	    /* If we move too far away from the initial step size then restart the Robbins-Munro
	       algorithm (Garthwaite, Fan & Sisson 2016) */
	    _theta0 = _theta;
	    _nstart = _nstart0 - _n;
	}

	_n++;
    }

    double StepAdapter::stepSize() const
    {
	return exp(_theta);
    }

    double StepAdapter::logitDeviation(double p) const
    {
	double logit_a = log(_a/(1 - _a));
	double logit_p = log(p/(1 - p));
	
	return logit_a - logit_p;
    }

} //namespace jags
