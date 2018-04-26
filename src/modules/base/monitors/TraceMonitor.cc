#include <config.h>
#include <graph/Node.h>

#include <algorithm>

#include "TraceMonitor.h"

using std::vector;
using std::string;

namespace jags {
namespace base {

    TraceMonitor::TraceMonitor(NodeArraySubset const &subset)
	: Monitor("trace", subset.nodes()), _subset(subset),
	  _values(subset.nchain())
    {
    }
    
    void TraceMonitor::update()
    {
	for (unsigned int ch = 0; ch < _values.size(); ++ch) {
	    vector<double> v = _subset.value(ch);
	    _values[ch].insert(_values[ch].end(), v.begin(), v.end());
	}
    }

    vector<double> const &TraceMonitor::value(unsigned int chain) const
    {
	return _values[chain];
    }

    vector<unsigned long> TraceMonitor::dim() const
    {
	return _subset.dim();
    }

    bool TraceMonitor::poolChains() const
    {
	return false;
    }

    bool TraceMonitor::poolIterations() const
    {
	return false;
    }

}}
