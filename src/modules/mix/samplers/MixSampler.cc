#include <config.h>
#include "MixSampler.h"
#include <graph/Graph.h>
#include <graph/StochasticNode.h>
#include <graph/GraphMarks.h>
#include <distribution/Distribution.h>
#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <sampler/GraphView.h>
#include <module/ModuleError.h>

#include <cmath>
#include <algorithm>

using std::vector;
using std::copy;
using std::logic_error;
using std::invalid_argument;
using std::log;
using std::exp;
using std::min;
using std::string;

#define N_REFRESH 100

namespace jags {

static void read_bounds(vector<StochasticNode*> const &snodes, 
			unsigned int chain,
			double *lower, double *upper, unsigned int length)
{

    double *lp = lower;
    double *up = upper;
    unsigned int node_length = 0;
    for (unsigned int j = 0; j < snodes.size(); ++j) {
	unsigned int length_j = snodes[j]->length();
	node_length += length_j;
	if (node_length > length) {
	    throwLogicError("Invalid length in read_bounds (MixSampler)");
	}
	else {
	    support(lp, up, length_j, snodes[j], chain);
	    lp += length_j;
	    up += length_j;
	}
    }
}

static vector<double> initialValue(GraphView const *gv, unsigned int chain)
{
    vector<double> ivalue(gv->length());
    gv->getValue(ivalue, chain);
    return ivalue;
}

namespace mix {

    MixSampler::MixSampler(GraphView const *gv, unsigned int chain,
			   unsigned int max_level, double delta, 
			   unsigned int nrep)
	: Metropolis(initialValue(gv, chain)),
	  _gv(gv), _chain(chain),
	  _max_level(max_level), 
	  _delta(delta),
	  _nrep(nrep),
	  _level(0),
	  _lower(0), _upper(0),
	  _n(0), 
	  _lstep(max_level+1), 
	  _pmean(max_level+1),
	  _temper(false)
    {
	if (_delta <= 0) {
	    throwLogicError("delta must be positive in MixSampler");
	}

	_lstep[0] = -5;
	_pmean[0] = 0;

	unsigned int N = _gv->length();
	_lower = new double[N];
	_upper = new double[N];
	read_bounds(_gv->nodes(), 0, _lower, _upper, N);
    }

    MixSampler::~MixSampler()
    {
	delete [] _lower;
	delete [] _upper;
    }

    //FIXME! This transformed random walk is useless without
    //Jacobian penalty!
    void MixSampler::updateStep(vector<double>  &x, double step, RNG *rng)
    {
	for (unsigned int i = 0; i < x.size(); ++i) {
	    bool bb = jags_finite(_lower[i]); //bounded below
	    bool ba = jags_finite(_upper[i]); //bounded above
	    double eps = rng->normal() * step;
	    if (bb && ba) {
		x[i] = log(x[i] - _lower[i]) - log(_upper[i] - x[i]);
		x[i] += eps;
		double w = 1.0 / (1 + exp(-x[i]));
		x[i] = (1 - w) * _lower[i] + w * _upper[i];
	    }
	    else if (bb) {
		x[i] = log(x[i] - _lower[i]); 
		x[i] += eps;
		x[i] = _lower[i] + exp(x[i]);
	    }
	    else if (ba) {
		x[i] = log(_upper[i] - x[i]);
		x[i] += eps;
		x[i] = _upper[i] - exp(x[i]);
	    }
	    else {
		x[i] = x[i] + eps;
	    }
	}
    }

    bool MixSampler::canSample(vector<StochasticNode *> const &snodes,
			       Graph const &graph)
    {
	if (snodes.empty())
	    return false;

	for (unsigned int i = 0; i < snodes.size(); ++i) {
	    if (!graph.contains(snodes[i])) {
		return false;
	    }
	    if (snodes[i]->isDiscreteValued()) {
		return false;
	    }
	    /* Check that the bounds (if any) of each node are fixed 
	       We rely on parameter transformation for bounded nodes which
	       must remain the same between iterations to guarantee
	       stationarity of the chain
	    */
	    if (!isSupportFixed(snodes[i])) {
		return false;
	    }
	    //Check that all nodes are of full rank
	    if (!snodes[i]->fullRank())
		return false;
	}
	return true;
    }

    void MixSampler::update(RNG *rng)
    {
	unsigned int length = _gv->length();
	vector<double> value1(length);
	getValue(value1);

	//We first do a non-tempered update
	
	_temper = false;
	double pmean = 0;
	for (unsigned int i = 0; i < _nrep; ++i) {
	    double lprob = -_gv->logFullConditional(_chain);
	    getValue(value1);
	    updateStep(value1, exp(_lstep[0]),rng);
	    setValue(value1);
	    lprob += _gv->logFullConditional(_chain);
	    double prob = min(exp(lprob), 1.0);
	    accept(rng, prob);
	    pmean += prob / _nrep;
	}
	_pmean[0] += pmean / N_REFRESH;


	//Now do a tempered update

	_temper = true;
    
	unsigned int nstep = 2 * _level;
	vector<double> pwr(nstep + 2);
	for (unsigned int t = 0; t <= _level; ++t) {
	    // NB Operator precedence here! Since t is an unsigned integer
	    // you get something very nasty if you try exp(-t * _delta)
	    pwr[t] = pwr[nstep + 1 - t] = exp(-(t * _delta));
	}


	//Save a copy of the current state(value, likelihood, prior)
	vector<double> value0(length);
	getValue(value0);
	double lprior0 = _gv->logPrior(_chain);
	double llik0 = _gv->logLikelihood(_chain);

	double log_global_prob = (pwr[1] - pwr[0]) * llik0;
	
	for (unsigned int t = 1; t <= nstep; ++t) {

	    unsigned int l = (t <= _level) ? t : (nstep + 1 - t);
	    double step = exp(_lstep[l]);

	    pmean = 0;
	    for (unsigned int i = 0; i < _nrep; ++i) {
		// Generate new proposal
		copy(value0.begin(), value0.end(), value1.begin());
		updateStep(value1, step, rng);
		setValue(value1);
	    
		// Calculate new prior and likelihood
		double lprior1 = _gv->logPrior(_chain);
		double llik1 = _gv->logLikelihood(_chain);
	    
		// Calculate acceptance probability for new proposal
		double lprob = (lprior1 - lprior0) + pwr[t] * (llik1 - llik0);
		double prob = exp(lprob);
		if (rng->uniform() <= prob) {
		    //Accept modified proposal
		    copy(value1.begin(), value1.end(), value0.begin());
		    lprior0 = lprior1;
		    llik0 = llik1;
		}
		pmean += min(prob, 1.0) / _nrep;
	    }
	    log_global_prob += (pwr[t+1] - pwr[t]) * llik0;
	    if (t <= _level) {
		_pmean[t] += pmean / N_REFRESH;
	    }
	}
	//Revert to last known good proposal
	setValue(value0); 
	accept(rng, exp(log_global_prob));
    }

    void MixSampler::rescale(double prob)
    {
	if (_temper) {

	    /* Rescale only after a tempered update */

	    _n++;
	    if (_n % N_REFRESH != 0) {
		return;
	    }

	    bool adapted = true;
	    for (unsigned t = 0; t <= _level; ++t) {
		if (_pmean[t] < 0.05) {
		    _lstep[t] -= log(1.50);
		}
		else if (_pmean[t] < 0.15) {
		    _lstep[t] -= log(1.15);
		}
		else if (_pmean[t] > 0.90) {
		    _lstep[t] += log(1.50);
		}
		else if (_pmean[t] > 0.35) {
		    _lstep[t] += log(1.15);
		}
		if (_pmean[t] < 0.15 || _pmean[t] > 0.35) {
		    adapted = false;
		}
		_pmean[t] = 0; //Reset
	    }
	    if (adapted && _level < _max_level) {
		_level++;
		_lstep[_level] = _lstep[_level - 1] + _delta/2;
		_pmean[_level] = 0;
	    }
	}
    }

    bool MixSampler::checkAdaptation() const
    {
	/* True if all up-phase transitions are within target range
	   for acceptance */
	return (_level == _max_level);
    }

    string MixSampler::name() const
    {
	return "MixSampler";
    }

    void MixSampler::getValue(vector<double> &x) const 
    {
	_gv->getValue(x, _chain);
    }

    void MixSampler::setValue(vector<double> const &x)
    {
	_gv->setValue(x, _chain);
    }

}}
