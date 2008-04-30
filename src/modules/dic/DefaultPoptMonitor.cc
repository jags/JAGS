#include "DefaultPoptMonitor.h"

#include <graph/StochasticNode.h>
#include <rng/RNG.h>

#include <stdexcept>
#include <cmath>

using std::vector;
using std::logic_error;
using std::string;
using std::exp;

namespace dic {
    
    DefaultPoptMonitor::DefaultPoptMonitor(StochasticNode const *snode,
					   unsigned int start, 
					   unsigned int thin, 
					   vector<RNG *> const &rngs, 
					   unsigned int nrep)
	: PoptMonitor(snode, start, thin), 
	  _repnode(snode->distribution(), snode->parents(), 
                   snode->lowerBound(), snode->upperBound()),
	  _rngs(rngs), _nrep(nrep)
	  
    {
    }


    void DefaultPoptMonitor::doUpdate()
    {
	unsigned int nchain = _repnode.nchain();
	unsigned int len = _repnode.length();

	vector<double> w(nchain); //weights
	double wsum = 0;
	for (unsigned int i = 0; i <  nchain; ++i) {
	    w[i] = exp(-_snode->logDensity(i));
	    wsum += w[i];
	    _weights[i] += w[i];
	}

	double pdsum = 0;
	for (unsigned int r = 0; r < _nrep; ++r) {
	    
	    for (unsigned int i = 0; i < nchain; ++i) {
		_repnode.randomSample(_rngs[i], i);
		double loglik = (wsum - w[i]) * _repnode.logDensity(i);
		
		for (unsigned int j = 0; j < nchain; ++j) {
		    if (j != i) {
			_repnode.setValue(_repnode.value(i), len, j);
			loglik -= w[j] * _repnode.logDensity(j);
		    }
		}
		pdsum += 2 * w[i] * loglik;
	    }
	}
	_values.push_back(pdsum /  _nrep);
    }

}
