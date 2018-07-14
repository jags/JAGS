#include <config.h>
#include "MNormal.h"

#include <matrix.h>
#include <DMNorm.h>
#include <lapack.h>

#include <graph/StochasticNode.h>
#include <sampler/SingletonGraphView.h>
#include <rng/RNG.h>

#include <cmath>
#include <algorithm>

#define N_REFRESH 100

using std::vector;
using std::copy;
using std::exp;
using std::sqrt;
using std::min;
using std::string;

namespace jags {
namespace bugs {

static vector<double> initValue(SingletonGraphView const *gv, 
				unsigned int chain)
{
    double const *x = gv->node()->value(chain);
    unsigned long N = gv->node()->length();
    vector<double> ivalue(N);
    copy(x, x + N, ivalue.begin());
    return ivalue;
}

MNormMetropolis::MNormMetropolis(SingletonGraphView const *gv, 
				 unsigned int chain)
    : Metropolis(initValue(gv, chain)),
      _gv(gv), _chain(chain), 
      _mean(nullptr), _var(nullptr), _prec(nullptr), 
      _n(0), _n_isotonic(0), _sump(0), _meanp(0), _lstep(0), _nstep(10), 
      _p_over_target(true)
{
    gv->checkFinite(chain); //Check validity of initial values
    
    unsigned int N = gv->length();

    _mean = new double[N];
    _var = new double[N * N];
    _prec = new double[N * N];
    
    // Give initial mean zero and identity variance matrix 
    for (unsigned int i = 0; i < N; ++i) {
	_mean[i] = 0;
	for (unsigned int j = 0; j < N; ++j) {
	    _prec[i + N * j] = _var[i + N * j] = (i == j) ? 1 : 0;
	}
    }
}

MNormMetropolis::~MNormMetropolis()
{
    delete [] _mean;
    delete [] _var;
    delete [] _prec;
}

void MNormMetropolis::update(RNG *rng)
{
    double logdensity = -_gv->logFullConditional(_chain);
    double step = exp(_lstep);

    double const *xold = _gv->node()->value(_chain);
    unsigned long N = _gv->length();

    vector<double> eps(N);
    DMNorm::randomsample(&eps[0], nullptr, _var, false, N, rng);
    vector<double> xnew(N);
    for (unsigned int i = 0; i < N; ++i) {
	xnew[i] = xold[i] + eps[i] * step;
    }

    setValue(xnew);
    logdensity += _gv->logFullConditional(_chain);
    accept(rng, exp(logdensity));
}

void MNormMetropolis::rescale(double p)
{
    ++_n;
    p = min(p, 1.0);
    _sump += p;

    if (_n % N_REFRESH == 0) {
	//Calculate the running mean acceptance rate 
	_meanp = _sump / N_REFRESH;    
	_sump = 0;
    }

    if (_n_isotonic == 0) {
	//Adjust scale of proposal distribution to get optimal acceptance
	//rate using a noisy gradient algorithm
	_lstep += (p - 0.234) / _nstep;
	if ((p > 0.234) != _p_over_target) {
	    _p_over_target = !_p_over_target;
	    ++_nstep;
	}
	/* 
	   Isotonic random walk. Use the identity matrix (scaled by
	   the _lstep parameter) as the precision of the proposal
	   distribution until the acceptance rate lies in an interval
	   around the optimum.
	*/
	if (_n % N_REFRESH == 0 && _meanp >= 0.15 && _meanp <= 0.35) {
	    _n_isotonic = _n;
	    _nstep = 100; //reset the step size as we adapt proposal
	}
    }
    else {
        //This give better adaptation in the orange tree example
	_lstep += (p - 0.234) / sqrt(static_cast<double>(_nstep));
        _nstep++;
	/* 
	   Adaptive random walk: The variance of the proposal
	   distribution is adapted to the empirical variance of the
	   posterior distribution.
	   
	   We use weighted moment estimators for the mean and variance
	   so that more recent iterations get greter weight. This is
	   important because the chain has not converged: the effect
	   of an initial transient must be minimized.

	   The weights are proportional to (_n - _n_isotonic) for the
	   mean, and the iterative formula is exact.  For the
	   variance, the weights are proportional to _n, and the
	   formula is asymptotically correct.

	   For small values of (_n - _n_isotonic), the variance
	   estimator is shrunk towards the prior (identity matrix) in
	   order to ensure a smooth transition from the isotonic
	   random walk.
	*/

	unsigned int N = _gv->length();
	double const *x = _gv->node()->value(_chain);
	for (unsigned int i = 0; i < N; ++i) {
	    _mean[i] += 2 * (x[i] - _mean[i]) / (_n - _n_isotonic + 1);
	}
	for (unsigned int i = 0; i < N; ++i) {
	    for (unsigned int j = 0; j < N; ++j) {
		_var[i + N * j] += 2 * ((x[i] - _mean[i]) * 
					    (x[j] - _mean[j]) -
					    _var[i + N * j]) / _n;
	    }
	}
    }
}

bool MNormMetropolis::checkAdaptation() const
{
    return (_n_isotonic > 0) && (_meanp >= 0.15) && (_meanp <= 0.35);
}

void MNormMetropolis::getValue(vector<double> &value) const
{
    double const *v = _gv->node()->value(_chain);
    copy(v, v + _gv->length(), value.begin());
}

void MNormMetropolis::setValue(vector<double> const &value)
{
    _gv->setValue(value, _chain);
}

}}
