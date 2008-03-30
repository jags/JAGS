#include "DefaultPDMonitor.h"

#include <graph/StochasticNode.h>
#include <rng/RNG.h>

using std::vector;

namespace dic {

    DefaultPDMonitor::DefaultPDMonitor(StochasticNode const *snode,
			 unsigned int start, unsigned int thin, 
			 vector<RNG *> const &rngs, unsigned int nrep)
	: PDMonitor(snode, start, thin), 
	  _repnode(snode->distribution(), snode->parents(), 
                   snode->lowerBound(), snode->upperBound()),
           _rngs(rngs), _nrep(nrep)
    
    {
    }

    void DefaultPDMonitor::doUpdate()
    {
	unsigned int nchain = _repnode.nchain();
	unsigned int len = _repnode.length();
	
	double pdsum = 0;
	for (unsigned int r = 0; r < _nrep; ++r) {
	    for (unsigned int i = 0; i < nchain; ++i) {
		_repnode.randomSample(_rngs[i], i);
		pdsum += _repnode.logDensity(i);

		double const *v = _repnode.value(i);
		for (unsigned int j = 0; j < nchain; ++j) {
		    if (j != i) {
			_repnode.setValue(v, len, j);
			pdsum -= _repnode.logDensity(j)/(nchain - 1);
		    }
		}
	    }
	}
	_values.push_back(pdsum/(_nrep * nchain));
    }
}
