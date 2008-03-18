#include "PDMonitor.h"

#include <graph/StochasticNode.h>
#include <rng/RNG.h>

#include <stdexcept>

using std::vector;
using std::logic_error;
using std::string;

namespace dic {

    PDMonitor::PDMonitor(StochasticNode const *snode,
			 unsigned int start, unsigned int thin, 
			 vector<RNG *> const &rngs, unsigned int nrep)
	: Monitor("pD", snode, start, thin), 
	  _repnode(snode->distribution(), snode->parents(), 
                   snode->lowerBound(), snode->upperBound()),
           _rngs(rngs), _nrep(nrep)
    
    {
	if (snode->nchain() < 2) {
	    throw logic_error("PDMonitor needs at least 2 chains");
	}
    }

    unsigned int PDMonitor::nchain() const
    {
	return 1;
    }

    vector<unsigned int> PDMonitor::dim() const
    {
	return vector<unsigned int> (1,niter());
    }
 
    vector<double> const &PDMonitor::value(unsigned int chain) const
    {
	return _values;
    }

    void PDMonitor::doUpdate()
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

    void PDMonitor::reserve(unsigned int niter)
    {
	unsigned int N = 1 + niter / thin();
	_values.reserve(_values.size() + N);
    }

    SArray PDMonitor::dump() const
    {
	SArray ans(dim());
	ans.setValue(_values);

        vector<string> names(1,string("iteration"));
	ans.setDimNames(names);
	return ans;
    }

}
