#include "PoptMonitor.h"

#include <graph/StochasticNode.h>
#include <rng/RNG.h>

#include <stdexcept>
#include <cmath>

using std::vector;
using std::logic_error;
using std::string;
using std::exp;

namespace dic {

    PoptMonitor::PoptMonitor(StochasticNode const *snode,
			     unsigned int start, unsigned int thin, 
			     vector<RNG *> const &rngs, unsigned int nrep)
	: Monitor("popt", snode, start, thin), _snode(snode),
	  _repnode(snode->distribution(), snode->parents(), 
                   snode->lowerBound(), snode->upperBound()),
	  _rngs(rngs), _nrep(nrep), _weights(snode->nchain(), 0)
    
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

    void PoptMonitor::doUpdate()
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
		
		double const *v = _repnode.value(i);
		double lik1 = 0;
		for (unsigned int j = 0; j < nchain; ++j) {
		    if (j != i) {
			_repnode.setValue(v, len, j);
			loglik -= w[j] * _repnode.logDensity(j);
		    }
		}
		pdsum += 2 * w[i] * loglik;
	    }
	}
	_values.push_back(pdsum /  _nrep);
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
