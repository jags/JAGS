#include "DevianceMonitor.h"

#include <graph/StochasticNode.h>

#include <stdexcept>
#include <algorithm>

using std::vector;
using std::logic_error;
using std::string;
using std::copy;

namespace dic {

    DevianceMonitor::DevianceMonitor(StochasticNode const *snode, 
				     unsigned int start, unsigned int thin)
	: Monitor("deviance", snode, start, thin), _values(snode->nchain()), 
	  _snode(snode)
    {
    }

    unsigned int DevianceMonitor::nchain() const
    {
	return _snode->nchain();
    }

    vector<unsigned int> DevianceMonitor::dim() const
    {
	vector<unsigned int> d(2);
	d[0] = niter();
	d[1] = nchain();
	return d;
    }
 
    vector<double> const &DevianceMonitor::value(unsigned int chain) const
    {
	return _values[chain];
    }

    void DevianceMonitor::doUpdate()
    {
	unsigned int nchain = _snode->nchain();
	for (unsigned int ch = 0; ch < nchain; ++ch) {
	    _values[ch].push_back(-2 * _snode->logDensity(ch));
	}
    }

    void DevianceMonitor::reserve(unsigned int niter)
    {
	unsigned int N = 1 + niter / thin();
	for (unsigned int ch = 0; ch < _values.size(); ++ch) {
	    _values[ch].reserve(_values[ch].size() + N);
	}
    }
    
    SArray DevianceMonitor::dump() const
    {
	unsigned int length = _values[0].size(); 
	unsigned int nchain = _values.size();

	vector<double> v(length * nchain);
	vector<double>::iterator p = v.begin();
	for (unsigned int ch = 0; ch < nchain; ++ch) {
	    p = copy(_values[ch].begin(), _values[ch].end(), p);
	}

	SArray ans(dim());
	ans.setValue(v);    
	
	vector<string> names(2);
	names[0] = "iteration";
	names[1] = "chain";
	ans.setDimNames(names);
	return(ans);
    }
}
