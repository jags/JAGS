#include <config.h>

#include "WAICMonitor.h"
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>

#include <algorithm>

using std::vector;
using std::string;
using std::copy;
using std::fill;

namespace jags {

    static vector<Node const *> toNodeVec(vector<StochasticNode const *> const &s)
    {
	vector<Node const *> ans(s.size());
	copy (s.begin(), s.end(), ans.begin());
	return ans;
    }

    namespace dic {

	WAICMonitor::WAICMonitor(vector<StochasticNode const *> const &snodes)
	    : Monitor("mean", toNodeVec(snodes)), _snodes(snodes),
	      _nchain(snodes[0]->nchain()),
	      _mlik(_nchain, vector<double>(snodes.size(), 0)),
	      _vlik(_nchain, vector<double>(snodes.size(), 0)),
	      _values(snodes.size(), 0),
	      _n(1)
	{
	}

	WAICMonitor::~WAICMonitor() 
	{
	}

	vector<unsigned int> WAICMonitor::dim() const
	{
	    return vector<unsigned int> (1, _snodes.size());
	}
 
	vector<double> const &WAICMonitor::value(unsigned int chain) const
	{
	    return _values;
	}

	bool WAICMonitor::poolChains() const
	{
	    return true;
	}

	bool WAICMonitor::poolIterations() const
	{
	    return true;
	}

	void WAICMonitor::update()
	{
	    fill(_values.begin(), _values.end(), 0);
	    for (unsigned int ch = 0; ch < _nchain; ++ch) {
		for (unsigned int k = 0; k < _snodes.size(); ++k) {
		    double delta = _snodes[k]->logDensity(ch, PDF_LIKELIHOOD) -
			_mlik[ch][k];

		    _mlik[ch][k] += delta/_n;
		    if (_n > 1) {
			_vlik[ch][k] *= static_cast<double>(_n - 2)/(_n - 1);
			_vlik[ch][k] += delta * delta / _n;
		    }
		    _values[k] += _vlik[ch][k] / _nchain;
		}
	    }
	    _n++;
	}

    }
}
