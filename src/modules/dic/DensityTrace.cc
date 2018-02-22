#include <config.h>
#include <graph/Node.h>
#include <distribution/Distribution.h>
// Required for PDFtype enum

#include "DensityTrace.h"

#include <cmath>

using std::vector;
using std::string;

namespace jags {
namespace dic {

    DensityTrace::DensityTrace(vector<Node const *> const &nodes,
		DensityType const density_type, string const &monitor_name)
	: Monitor(monitor_name, nodes), _nodes(nodes), _density_type(density_type), 
		_nchain(nodes[0]->nchain()), _values(nodes[0]->nchain())
    {
		if( _density_type != DENSITY && _density_type != LOGDENSITY && _density_type != DEVIANCE ) {
			throw std::logic_error("Unimplemented DensityType in NodeDensityMonitorFactory");
		}
    }

    void DensityTrace::update()
    {
		for (unsigned int ch = 0; ch < _nchain; ++ch) {
		    for (unsigned int i = 0; i < _nodes.size(); ++i) {
				double newval = _nodes[i]->logDensity(ch, PDF_FULL);
				if( _density_type == DENSITY ) {
					newval = std::exp(newval);
				}
				else if ( _density_type == DEVIANCE ) {
					newval = -2.0 * newval;
				}
			    _values[ch].push_back(newval);
		    }
		}
    }
	
    vector<double> const &DensityTrace::value(unsigned int chain) const
    {
	return _values[chain];
    }

    vector<unsigned int> DensityTrace::dim() const
    {
	return vector<unsigned int>(1, _nodes.size());
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
