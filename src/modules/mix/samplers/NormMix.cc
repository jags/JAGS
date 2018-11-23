#include <config.h>
#include "NormMix.h"
#include "DirichletInfo.h"
#include <rng/RNG.h>
#include <util/nainf.h>
#include <sampler/GraphView.h>
#include <distribution/Distribution.h>
#include <module/ModuleError.h>

#include <cmath>

class Graph;

using std::vector;
using std::log;
using std::exp;
using std::string;
using std::pair;
using std::set;

namespace jags {

static bool isDirch(StochasticNode const *snode)
{
    return snode->distribution()->name() == "ddirch";
}

static vector<double> initialValue(GraphView const *gv, unsigned int chain)
{
    vector<double> ivalue(gv->length());
    gv->getValue(ivalue, chain);
    return ivalue;
}

static void reflect (double &x, double LLIM=-30, double ULIM=30)
{
    //boundary reflection on log scale to to prevent overflow/underflow
    //when taking exponentials.

    while (x < LLIM || x > ULIM) {
	if (x < LLIM) {
	    x = 2*LLIM - x;
	}
	if (x > ULIM) {
	    x = 2*ULIM - x;
	}
    }
}

namespace mix {

    NormMix::NormMix(GraphView const *gv, unsigned int chain,
		     unsigned int nlevel, double max_temp, unsigned int nrep)
	: TemperedMetropolis(initialValue(gv, chain), nlevel, max_temp, nrep),
	  _gv(gv), _chain(chain)
    {
	unsigned int N = gv->length();
	_lower = new double[N];
	_upper = new double[N];

	double *lp = _lower;
	double *up = _upper;
	unsigned long start = 0;
        vector<StochasticNode*> const &snodes = _gv->nodes();
	for (unsigned long j = 0; j < snodes.size(); ++j) {
	    unsigned long length_j = snodes[j]->length();
	    if (isDirch(snodes[j])) {
		//Special rule for Dirichlet distribution to enforce
		//log transformation
		for (unsigned long k = 0; k < length_j; ++k) {
		    lp[k] = 0;
		    up[k] = JAGS_POSINF;
		}
		_di.push_back(new DirichletInfo(snodes[j], start, chain));
	    }
	    else {
		snodes[j]->support(lp, up, length_j, chain);
	    }
	    lp += length_j;
	    up += length_j;
	    start += length_j;
	    /*
	    if (start > N) {
		throwLogicError("Invalid length in NormMix constructor");
	    }
	    */
	}
    }

    NormMix::~NormMix()
    {
	delete [] _lower;
	delete [] _upper;
	for (unsigned int i = 0; i < _di.size(); ++i) {
	    delete _di[i];
	}
    }

    void NormMix::step(vector<double>  &x, double step, RNG *rng) const
    {
	for (unsigned int i = 0; i < x.size(); ++i) {
	    bool bb = jags_finite(_lower[i]); //bounded below
	    bool ba = jags_finite(_upper[i]); //bounded above
	    double eps = rng->normal() * step;
		
	    if (bb && ba) {
		x[i] = log(x[i] - _lower[i]) - log(_upper[i] - x[i]);
		x[i] += eps;
		reflect(x[i]);
		double w = 1.0 / (1 + exp(-x[i]));
		x[i] = (1 - w) * _lower[i] + w * _upper[i];
	    }
	    else if (bb) {
		x[i] = log(x[i] - _lower[i]);
		x[i] += eps;
		reflect(x[i]);
		x[i] = _lower[i] + exp(x[i]);
	    }
	    else if (ba) {
		x[i] = log(_upper[i] - x[i]);
		x[i] += eps;
		reflect(x[i]);
		x[i] = _upper[i] - exp(x[i]);
	    }
	    else {
		x[i] = x[i] + eps;
	    }
	}
    }

    double NormMix::logJacobian(vector<double> const &x) const
    {
	double lj = 0;
	for (unsigned int i = 0; i < x.size(); ++i) {
	    bool bb = jags_finite(_lower[i]); //bounded below
	    bool ba = jags_finite(_upper[i]); //bounded above
	    if (bb && ba) {
		lj += log(x[i] - _lower[i]) + log(_upper[i] - x[i]);
	    }
	    else if (bb) {
		lj += log(x[i] - _lower[i]); 
	    }
	    else if (ba) {
		lj += log(_upper[i] - x[i]);
	    }
	}
	return lj;
    }

    bool NormMix::canSample(vector<StochasticNode *> const &snodes)
    {
	if (snodes.empty())
	    return false;

	for (unsigned int i = 0; i < snodes.size(); ++i) {
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
	    if (isDirch(snodes[i])) {
		/* We can't handle structural zeros in the Dirichlet
		   distribution. This excludes Dirichlet nodes with
		   non-fixed parents
		*/
		if (!snodes[i]->parents()[0]->isFixed()) {
		    return false;
		}
		double const *par = snodes[i]->parents()[0]->value(0);
		unsigned long plen = snodes[i]->parents()[0]->length();
		for (unsigned long j = 0; j < plen; ++j) {
		    if (par[j] == 0) {
			return false;
		    }
		}
	    }
	    else if (!snodes[i]->fullRank()) {
		return false;
	    }
	}
	return true;
    }

    void NormMix::getValue(vector<double> &x) const 
    {
	_gv->getValue(x, _chain);
	for (unsigned int i = 0; i < _di.size(); ++i) {
	    for (unsigned int j = _di[i]->start; j < _di[i]->end; ++j) {
		x[j] *= _di[i]->sum;
	    }
	}
    }

    void NormMix::setValue(vector<double> const &x)
    {
	if (_di.empty()) {
	    _gv->setValue(x, _chain);
	}
	else {
	    //Rescale Dirichlet distributions
	    for (unsigned int i = 0; i < _di.size(); ++i) {
		_di[i]->sum = 0;
		for (unsigned int j = _di[i]->start; j < _di[i]->end; ++j) {
		    _di[i]->sum += x[j];
		}
	    }
	    vector<double> v = x;
	    for (unsigned int i = 0; i < _di.size(); ++i) {
		for (unsigned int j = _di[i]->start; j < _di[i]->end; ++j) {
		    v[j] /= _di[i]->sum;
		}
	    }
	    _gv->setValue(v, _chain);
	}
    }

    double NormMix::logPrior() const
    {
	double lp = _gv->logPrior(_chain);
	for (unsigned int i = 0; i < _di.size(); ++i) {
	    //Penalty on unscaled Dirichlet coordinates
	    //Gamma penalty
	    lp += _di[i]->gammaPenalty();
	    /*
	      double logS = log(_di[i]->sum);
	      lp -= 5 * logS * logS;
	    */
	}
	return lp;
    }
    
    double NormMix::logLikelihood() const
    {
        return _gv->logLikelihood(_chain);
    }

}}
