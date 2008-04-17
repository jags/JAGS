#include "KLPDMonitor.h"

#include <graph/StochasticNode.h>
#include <rng/RNG.h>

using std::vector;

namespace dic {

    KLPDMonitor::KLPDMonitor(StochasticNode const *snode,
			     unsigned int start, unsigned int thin,
			     KL const *kl)
	: PDMonitor(snode, start, thin), _kl(kl)
    {
    }

    void KLPDMonitor::doUpdate()
    {
	unsigned int nchain = _snode->nchain();
	
	double pdsum = 0;
	for (unsigned int i = 1; i < nchain; ++i) {
	    for (unsigned int j = 0; j < i; ++j) {
		pdsum += _kl->divergence(_snode->parameters(i),
					 _snode->parameters(j));
	    }
	}

	_values.push_back(pdsum/(nchain * (nchain - 1)));
    }
}
