#include <config.h>
#include <graph/Node.h>
#include <util/nainf.h>

#include <algorithm>

#include "PoolVarianceMonitor.h"

using std::vector;
using std::string;

namespace jags {
namespace base {

    PoolVarianceMonitor::PoolVarianceMonitor(NodeArraySubset const &subset)
	: Monitor("poolvariance", subset.nodes()), _subset(subset),
	  _means(subset.length()),
	  _mms(subset.length()),
	  _variances(subset.length()),
	  _n(0)
    {
    }
    
    void PoolVarianceMonitor::update()
    {

		for (unsigned int ch = 0; ch < _subset.nchain(); ++ch) {
		
			// Each chain counts as an iteration:
			_n++;
		
		    vector<double> value = _subset.value(ch);
		    for (unsigned int i = 0; i < value.size(); ++i) {
				if (value[i] == JAGS_NA) {
				    _means[i] = JAGS_NA;
					_mms[i] = JAGS_NA;
					_variances[i] = JAGS_NA;
				}
				else {
					double delta = value[i] - _means[i];
					_means[i] += delta / _n;
					_mms[i] += delta * (value[i] - _means[i]);
				}
			}
		}
		
		// Variance itself only needs to be calculated once per iteration:
		for (unsigned int i = 0; i < _variances.size(); ++i) {
		    _variances[i] = _mms[i] / static_cast<double>(_n - 1);
		}
		
    }

    vector<double> const &PoolVarianceMonitor::value(unsigned int) const
    {
	return _variances;
    }
    
    vector<unsigned long> PoolVarianceMonitor::dim() const
    {
	return _subset.dim();
    }

    bool PoolVarianceMonitor::poolChains() const
    {
	return true;
    }

    bool PoolVarianceMonitor::poolIterations() const
    {
	return true;
    }
	
}}
