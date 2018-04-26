#include <config.h>
#include <graph/Node.h>
#include <util/nainf.h>

#include <algorithm>

#include "PoolMeanMonitor.h"

using std::vector;
using std::string;

namespace jags {
namespace base {

    PoolMeanMonitor::PoolMeanMonitor(NodeArraySubset const &subset)
	: Monitor("poolmean", subset.nodes()), _subset(subset),
	  _values(subset.length()),
	  _n(0)
    {
	
    }
    
    void PoolMeanMonitor::update()
    {

		for (unsigned int ch = 0; ch < _subset.nchain(); ++ch) {

			// Each chain counts as an iteration:
			_n++;

		    vector<double> value = _subset.value(ch);
		    for (unsigned int i = 0; i < value.size(); ++i) {
				if (value[i] == JAGS_NA) {
				    _values[i] = JAGS_NA;
				}
				else {
				    _values[i] -= (_values[i] - value[i])/_n;
				}
			}
		}
    }

    vector<double> const &PoolMeanMonitor::value(unsigned int chain) const
    {
	return _values;
    }

    vector<unsigned long> PoolMeanMonitor::dim() const
    {
	return _subset.dim();
    }

    bool PoolMeanMonitor::poolChains() const
    {
	return true;
    }

    bool PoolMeanMonitor::poolIterations() const
    {
	return true;
    }

}}
