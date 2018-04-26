#include "DevianceTrace.h"

#include <graph/StochasticNode.h>

#include <algorithm>

using std::vector;
using std::string;
using std::copy;

namespace jags {

static vector<Node const *> toNode(vector<StochasticNode const *> const &nodes)
{
    vector<Node const *> ans(nodes.size());
    for (unsigned int i = 0; i < nodes.size(); ++i) {
	ans[i] = nodes[i];
    }
    return ans;
}

namespace dic {

    DevianceTrace::DevianceTrace(vector<StochasticNode const *> const &
				     snodes)
	: Monitor("trace", toNode(snodes)), _values(snodes[0]->nchain()), 
	  _snodes(snodes)
    {
    }

    vector<unsigned long> DevianceTrace::dim() const
    {
	return vector<unsigned long>(1,1);
    }
 
    vector<double> const &DevianceTrace::value(unsigned int chain) const
    {
	return _values[chain];
    }

    void DevianceTrace::update()
    {
	unsigned int nchain = _snodes[0]->nchain();
	for (unsigned int ch = 0; ch < nchain; ++ch) {
	    double loglik = 0;
	    for (unsigned long i = 0; i < _snodes.size(); ++i) {
		loglik += _snodes[i]->logDensity(ch, PDF_FULL);
	    }
	    _values[ch].push_back(-2 * loglik);
	}
    }

    bool DevianceTrace::poolChains() const
    {
	return false;
    }

    bool DevianceTrace::poolIterations() const
    {
	return false;
    }

}}
