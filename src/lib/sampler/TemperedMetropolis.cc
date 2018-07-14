#include <config.h>
#include <sampler/TemperedMetropolis.h>
#include <rng/RNG.h>

#include <cmath>
#include <stdexcept>

using std::vector;
using std::invalid_argument;
using std::log;
using std::exp;
using std::fabs;

//Minimum number of iterations before we can go to the next level
#define MIN_STEP 50

static vector<double> makePower(unsigned int max_level, double max_temp)
{
    vector<double> pwr(max_level + 1);
    double delta = log(max_temp) / max_level;
    for (unsigned int t = 0; t <= max_level; ++t) {
	pwr[t] = exp(-t * delta);
    }
    return pwr;
}

namespace jags {

TemperedMetropolis::TemperedMetropolis(vector<double> const &value,
				       unsigned int max_level, double max_temp, 
				       unsigned int nrep)
    : Metropolis(value),
      _max_level(max_level), 
      _nrep(nrep),
      _pwr(makePower(max_level, max_temp)),
      _t(0), _tmax(1), _step_adapter(0),
      _pmean(0), _niter(2)
{
    if (max_temp <= 1) {
	throw invalid_argument("Invalid max_temp in TemperedMetropolis ");
    }
    if (max_level <= 0) {
	throw invalid_argument("Invalid max_level in TemperedMetropolis");
    }
    
    _step_adapter.push_back(nullptr);
    StepAdapter *adapter = new StepAdapter(0.1);
    _step_adapter.push_back(adapter);
}

TemperedMetropolis::~TemperedMetropolis()
{
    for (unsigned int i = 1; i < _step_adapter.size(); ++i) {
	delete _step_adapter[i];
    }
}

void TemperedMetropolis::temperedUpdate(RNG *rng, 
					double &lprior0,
					double &llik0, 
					vector<double> &value0)
{
    vector<double> x = value0;
    for (unsigned int r = 0; r < _nrep; ++r) {
	step(x, _step_adapter[_t]->stepSize(), rng);
	setValue(x);
	double lprior1 = logPrior() + logJacobian(x);
	double llik1 = logLikelihood();
	double lprob = (lprior1 - lprior0) +  _pwr[_t] * (llik1 - llik0);
	if (accept(rng, exp(lprob))) {
	    lprior0 = lprior1;
	    llik0 = llik1;
	    value0 = x;
	}
	else {
	    x = value0;
	}
    }
}

void TemperedMetropolis::update(RNG *rng)
{
    //Save the current state
    vector<double> last_value(length());
    getValue(last_value);

    //Make copies of the current state and log density
    //These are modified in place by temperedUpdate
    double log_prior = logPrior() + logJacobian(last_value);
    double log_lik = logLikelihood();
    vector<double> current_value = last_value;

    double log_global_prob = 0;
    for (_t = 1U; _t <= _tmax; _t++) {
	log_global_prob += (_pwr[_t] - _pwr[_t-1]) * log_lik;
	temperedUpdate(rng, log_prior, log_lik, current_value);
    }
    for (_t = _tmax; _t >= 1; _t--) {
	temperedUpdate(rng, log_prior, log_lik, current_value);
	log_global_prob -= (_pwr[_t] - _pwr[_t-1])  * log_lik;
    }

    //Global Metropolis-Hastings acceptance step
    if (!accept(rng, exp(log_global_prob))) {
	//Force return to last good value
	setValue(last_value);
	accept(rng, 1.0);
    }
}

void TemperedMetropolis::rescale(double p)
{
    if (_t == 0)
	return; //No adaptation for global acceptance step

    _step_adapter[_t]->rescale(p);

    if (_t == _tmax && _tmax < _max_level) {

	// We keep a weighted mean estimate of the mean acceptance probability
	//  with the weights in favour of more recent iterations
	_pmean += 2 * (p - _pmean) / _niter++;

	// If mean acceptance probability is in target range, increase
	// temperature to the next level. We enforce a delay of
	// MIN_STEP iterations first.
	double delta = _step_adapter[_t]->logitDeviation(_pmean);
	if (_niter > MIN_STEP + 2 && fabs(delta) < 0.25) {
	    _niter = 2;
	    _pmean = 0;
	    _tmax++;
	    double step = _step_adapter.back()->stepSize();
	    StepAdapter *adapter = new StepAdapter(step);
	    _step_adapter.push_back(adapter);
	}
    }
}


bool TemperedMetropolis::checkAdaptation() const
{
    return (_tmax == _max_level);
}

void TemperedMetropolis::step(vector<double> &value, double s, RNG *rng) const
{
    for (unsigned int i = 0; i < value.size(); ++i) {
	value[i] += rng->normal() * s;
    }
}

double TemperedMetropolis::logJacobian(vector<double> const &) const
{
    return 0;
}

} //namespace jags
