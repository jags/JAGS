#include "PoptMonitor.h"

#include <graph/StochasticNode.h>

#include <stdexcept>
#include <cmath>

using std::vector;
using std::logic_error;
using std::string;
using std::exp;

namespace dic {

    PoptMonitor::PoptMonitor(StochasticNode const *snode,
			     unsigned int start, unsigned int thin) 
	: Monitor("popt", snode, start, thin), _snode(snode),
	  _weights(snode->nchain(), 0)
    
    {
	if (snode->nchain() < 2) {
	    throw logic_error("PoptMonitor needs at least 2 chains");
	}
    }

    unsigned int PoptMonitor::nchain() const
    {
	return 1;
    }

    vector<unsigned int> PoptMonitor::dim() const
    {
	return vector<unsigned int> (1,niter());
    }
 
    vector<double> const &PoptMonitor::value(unsigned int chain) const
    {
	return _values;
    }

    void PoptMonitor::reserve(unsigned int niter)
    {
	unsigned int N = 1 + niter / thin();
	_values.reserve(_values.size() + N);
    }

    SArray PoptMonitor::dump() const
    {
	SArray ans(dim());
	vector<double> scaled_values(_values);

	double wsum = 0;
	for (unsigned int i = 0; i < _weights.size(); ++i) {
	   for (unsigned int j = 0; j < _weights.size(); ++j) {
              if (j != i) {
                 wsum += _weights[i] * _weights[j];
              }
           }
	}

        //double scale = niter() / (wsum2 - wsum * wsum);
        double scale = niter() * niter() / wsum;
	for (unsigned int i = 0; i < _values.size(); ++i) {
	    scaled_values[i] *= scale;
	}
	ans.setValue(scaled_values);

        vector<string> names(1,string("iteration"));
	ans.setDimNames(names);
	return ans;
    }

}
