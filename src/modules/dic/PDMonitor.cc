#include <config.h>

#include "PDMonitor.h"
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>

#include <algorithm>

using std::vector;
using std::string;
using std::copy;

namespace jags {

static vector<Node const *> toNodeVec(vector<StochasticNode const *> const &s)
{
    vector<Node const *> ans(s.size());
    copy (s.begin(), s.end(), ans.begin());
    return ans;
}

namespace dic {

    PDMonitor::PDMonitor(vector<StochasticNode const *> const &snodes,
			 vector<RNG *> const &rngs,
			 unsigned int nrep, double scale)
	: Monitor("mean", toNodeVec(snodes)), _snodes(snodes), _rngs(rngs),
	  _nrep(nrep),
	  _values(snodes.size(), 0),  _weights(snodes.size(), 0),
	  _scale(scale), _nchain(rngs.size())
    {
	if (_nchain < 2) {
	    throwLogicError("PDMonitor needs at least 2 chains");
	}
    }

    PDMonitor::~PDMonitor() 
    {
    }

    vector<unsigned long> PDMonitor::dim() const
    {
	return vector<unsigned long> (1,_values.size());
    }
 
    vector<double> const &PDMonitor::value(unsigned int chain) const
    {
	return _values;
    }

    bool PDMonitor::poolChains() const
    {
	return true;
    }

    bool PDMonitor::poolIterations() const
    {
	return true;
    }

    void PDMonitor::update()
    {
	vector<double> w(_nchain);
	for (unsigned int k = 0; k < _values.size(); ++k) {
	    
	    double pdsum = 0;
	    double wsum = 0;
	    for (unsigned int i = 0; i < _nchain; ++i) {
		w[i] = weight(_snodes[k], i);
		for (unsigned int j = 0; j < i; ++j) {
		    pdsum += w[i] * w[j] * (
			_snodes[k]->KL(i, j, _rngs[i], _nrep) +
			_snodes[k]->KL(j, i, _rngs[j], _nrep));
		    wsum += w[i] * w[j];
		}
	    }
	    pdsum /= wsum;
	    pdsum *= _scale/2;

	    _weights[k] += wsum;
	    _values[k] += wsum * (pdsum - _values[k])/_weights[k];
	}
    }


    double PDMonitor::weight(StochasticNode const *snode, unsigned int ch) const
    {
	return 1;
    }

}}
