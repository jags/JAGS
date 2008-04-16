#include "KLPoptMonitor.h"

#include <graph/StochasticNode.h>
#include <rng/RNG.h>

#include <stdexcept>
#include <cmath>

using std::vector;
using std::logic_error;
using std::string;
using std::exp;

namespace dic {

    KLPoptMonitor::KLPoptMonitor(StochasticNode const *snode,
			     unsigned int start, unsigned int thin, 
			     vector<RNG *> const &rngs, unsigned int nrep)
	: PoptMonitor(snode, start, thin)
    
    {
	unsigned int nchain = snode->nchain();
	_par.reserve(nchain);
	for (unsigned int i = 0; i < nchain; ++i) {
	    _par.push_back(snode->parameters(i));
	}
    }

    KLPoptMonitor::~KLPoptMonitor()
    {
	delete _kl;
    }

    void KLPoptMonitor::doUpdate()
    {
	unsigned int nchain = _par.size();

	vector<double> w(nchain); //weights
	double wsum = 0;
	for (unsigned int i = 0; i <  nchain; ++i) {
	    w[i] = exp(-_snode->logDensity(i));
	    wsum += w[i];
	    _weights[i] += w[i];
	}

	double pdsum = 0;
	for (unsigned int i = 0; i < nchain; ++i) {
	    double pdi = 0;
	    for (unsigned int j = 0; j < nchain; ++j) {
		if (j != i) {
		    pdsum += w[j] * _kl->divergence(_par[i], _par[j]);
		}
	    }
	    pdsum += 2 * w[i] * pdi;
	}

	_values.push_back(pdsum);
    }
}
