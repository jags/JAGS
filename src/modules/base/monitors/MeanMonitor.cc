#include <config.h>
#include <graph/Node.h>
#include <util/nainf.h>

#include <algorithm>

#include "MeanMonitor.h"

using std::vector;
using std::string;

namespace jags {
namespace base {

    MeanMonitor::MeanMonitor(NodeArraySubset const &subset)
	: Monitor("mean", subset.nodes()), _subset(subset),
	  _values(subset.nchain(), vector<double>(subset.length())),
	  _n(0)
    {
	
    }
    
    void MeanMonitor::update()
    {
	_n++;
	for (unsigned int ch = 0; ch < _values.size(); ++ch) {
	    vector<double> value = _subset.value(ch);
	    vector<double> &rmean  = _values[ch];
	    for (unsigned int i = 0; i < value.size(); ++i) {
		if (value[i] == JAGS_NA) {
		    rmean[i] = JAGS_NA;
		}
		else {
		    rmean[i] -= (rmean[i] - value[i])/_n;
		}
	    }
	}
    }

    vector<double> const &MeanMonitor::value(unsigned int chain) const
    {
	return _values[chain];
    }

    vector<unsigned long> MeanMonitor::dim() const
    {
	return _subset.dim();
    }

    bool MeanMonitor::poolChains() const
    {
	return false;
    }

    bool MeanMonitor::poolIterations() const
    {
	return true;
    }

}}
