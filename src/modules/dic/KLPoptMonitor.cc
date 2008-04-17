#include "KLPoptMonitor.h"

#include <graph/StochasticNode.h>
#include <rng/RNG.h>

#include <cmath>

using std::vector;
using std::string;
using std::exp;

namespace dic {

    KLPoptMonitor::KLPoptMonitor(StochasticNode const *snode,
				 unsigned int start, unsigned int thin, 
				 KL const *kl)
	: PoptMonitor(snode, start, thin), _kl(kl)
    
    {
    }

    void KLPoptMonitor::doUpdate()
    {
	unsigned int nchain = _snode->nchain();

	vector<double> w(nchain); //weights
	double wsum = 0;
	for (unsigned int i = 0; i <  nchain; ++i) {
	    w[i] = exp(-_snode->logDensity(i));
	    wsum += w[i];
	    _weights[i] += w[i];
	}

	double pdsum = 0;
	for (unsigned int i = 1; i < nchain; ++i) {
	    for (unsigned int j = 0; j < i; ++j) {
		pdsum += 2*w[i]*w[j] * _kl->divergence(_snode->parameters(i),
						       _snode->parameters(j));
	    }
	}

	_values.push_back(pdsum);
    }
}
