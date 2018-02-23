#include <config.h>
#include <graph/Node.h>
#include <distribution/Distribution.h>
// Required for PDFtype enum
#include <util/nainf.h>

#include "DensityTrace.h"

#include <cmath>

using std::vector;
using std::string;

namespace jags {
namespace dic {

    DensityTrace::DensityTrace(vector<Node const *> const &nodes, vector<unsigned int> dim,
		DensityType const density_type, string const &monitor_name)
	: Monitor(monitor_name, nodes), _nodes(nodes), _density_type(density_type), 
		_nchain(nodes[0]->nchain()), _values(nodes[0]->nchain()), _dim(dim)
    {
		if( _density_type != DENSITY && _density_type != LOGDENSITY && _density_type != DEVIANCE ) {
			throw std::logic_error("Unimplemented DensityType in DensityTrace");
		}
    }

    void DensityTrace::update()
    {
		for (unsigned int ch = 0; ch < _nchain; ++ch) {
		    for (unsigned int i = 0; i < _nodes.size(); ++i) {
				double newval = _nodes[i]->logDensity(ch, PDF_FULL);
				if (newval == JAGS_NA) {
				    // Don't try and convert NA to density or deviance
				}else if( _density_type == DENSITY ) {
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
	return _dim;
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
