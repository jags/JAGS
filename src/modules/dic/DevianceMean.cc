#include "DevianceMean.h"

#include <graph/StochasticNode.h>

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

    DevianceMean::DevianceMean(vector<StochasticNode const *> const &s)
	: Monitor("mean", toNodeVec(s)), _values(s.size(),0), _snodes(s), _n(0)
    {
    }

    vector<unsigned long> DevianceMean::dim() const
    {
	return vector<unsigned long>(1, _values.size());
    }

    bool DevianceMean::poolChains() const
    {
	return true;
    }

    bool DevianceMean::poolIterations() const
    {
	return true;
    }

    vector<double> const &DevianceMean::value(unsigned int chain) const
    {
	return _values;
    }

    void DevianceMean::update()
    {
	_n++;
	for (unsigned long i = 0; i < _snodes.size(); ++i) {
	    double loglik = 0;
	    unsigned int nchain = _snodes[i]->nchain();
	    for (unsigned int ch = 0; ch < nchain; ++ch) {
		loglik += _snodes[i]->logDensity(ch, PDF_FULL) / nchain;
	    }
	    _values[i] += (-2*loglik - _values[i])/_n;
	}
    }

}}
