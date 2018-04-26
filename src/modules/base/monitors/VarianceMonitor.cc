#include <config.h>
#include <graph/Node.h>
#include <util/nainf.h>

#include <algorithm>

#include "VarianceMonitor.h"

using std::vector;
using std::string;

namespace jags {
namespace base {

    VarianceMonitor::VarianceMonitor(NodeArraySubset const &subset)
	: Monitor("variance", subset.nodes()), _subset(subset),
	  _means(subset.nchain(), vector<double>(subset.length())),
	  _mms(subset.nchain(), vector<double>(subset.length())),
	  _variances(subset.nchain(), vector<double>(subset.length())),
	  _n(0)
    {
    }
    
    void VarianceMonitor::update()
    {
	_n++;
	for (unsigned int ch = 0; ch < _means.size(); ++ch) {
	    vector<double> value = _subset.value(ch);
	    vector<double> &rmean  = _means[ch];
	    vector<double> &rmm  = _mms[ch];
		vector<double> &rvar  = _variances[ch];		
	    for (unsigned int i = 0; i < value.size(); ++i) {
		if (value[i] == JAGS_NA) {
		    rmean[i] = JAGS_NA;
			rmm[i] = JAGS_NA;
			rvar[i] = JAGS_NA;
		}
		else {
			double delta = value[i] - rmean[i];
			rmean[i] += delta / _n;
			rmm[i] += delta * (value[i] - rmean[i]);
			
			/* 
			[MD - 2016-03-31]
			The variance would be more efficiently calculated on-demand,
			but the functions are all marked as const except this one.
			It is possible to get around this restriction, but the value
			function is const for a good reason, so hacking around this
			seems to me to be a bad idea.
			I do have an alternative implementation that calculates
			variance on demand, which we could consider using for the
			sake of efficiency...  
			*/
			rvar[i] = rmm[i] / static_cast<double>(_n - 1);
			
		}
	}
	}
    }

    vector<double> const &VarianceMonitor::value(unsigned int chain) const
    {
		return _variances[chain];
	}
	
    vector<unsigned long> VarianceMonitor::dim() const
    {
	return _subset.dim();
    }

    bool VarianceMonitor::poolChains() const
    {
	return false;
    }

    bool VarianceMonitor::poolIterations() const
    {
	return true;
    }
	
}}
