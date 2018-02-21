#include <config.h>
#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>
// Required for PDFtype enum

#include <algorithm>

#include "LogDensityMonitor.h"

using std::vector;
using std::string;

namespace jags {
namespace base {

    LogDensityMonitor::LogDensityMonitor(NodeArraySubset const &subset)
	: Monitor("logdensity", subset.nodes()), _subset(subset),
	  _values(subset.nchain())
    {
    }
    
    void LogDensityMonitor::update()
    {
	for (unsigned int ch = 0; ch < _values.size(); ++ch) {
		vector<double> v = _subset.logDensity(ch, PDF_FULL);
//	    vector<double> v = _subset.value(ch);
	    _values[ch].insert(_values[ch].end(), v.begin(), v.end());
	}
    }

    vector<double> const &LogDensityMonitor::value(unsigned int chain) const
    {
	return _values[chain];
    }

    vector<unsigned int> LogDensityMonitor::dim() const
    {
	return _subset.dim();
    }

    bool LogDensityMonitor::poolChains() const
    {
	return false;
    }

    bool LogDensityMonitor::poolIterations() const
    {
	return false;
    }

}}
