#include <config.h>
#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>
// Required for PDFtype enum

#include "DensityTrace.h"

using std::vector;
using std::string;

namespace jags {
namespace dic {

    DensityTrace::DensityTrace(NodeArraySubset const &subset, DensityType const density_type, string const &monitor_name)
	: Monitor(monitor_name, subset.nodes()), _subset(subset),
	  _values(subset.nchain()), _density_type(density_type)
    {
    }

    void DensityTrace::update()
    {
	for (unsigned int ch = 0; ch < _values.size(); ++ch) {
		vector<double> v = _subset.logDensity(ch, PDF_FULL, _density_type);
	    _values[ch].insert(_values[ch].end(), v.begin(), v.end());
	}
    }
	
    vector<double> const &DensityTrace::value(unsigned int chain) const
    {
	return _values[chain];
    }

    vector<unsigned int> DensityTrace::dim() const
    {
	return _subset.dim();
    }

    bool DensityTrace::poolChains() const
    {
	return false;
    }

    bool DensityTrace::poolIterations() const
    {
	return false;
    }

}}
