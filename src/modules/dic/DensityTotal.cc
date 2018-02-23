#include <config.h>
#include <graph/Node.h>
#include <distribution/Distribution.h>
// Required for PDFtype enum
#include <util/nainf.h>

#include "DensityTotal.h"

#include <cmath>

using std::vector;
using std::string;

namespace jags {
namespace dic {

    DensityTotal::DensityTotal(vector<Node const *> const &nodes, vector<unsigned int> dim,
		DensityType const density_type, string const &monitor_name)
	: Monitor(monitor_name, nodes), _nodes(nodes), _density_type(density_type), 
		_nchain(nodes[0]->nchain()), _values(nodes[0]->nchain())
    {
		if( _density_type != DENSITY && _density_type != LOGDENSITY && _density_type != DEVIANCE ) {
			throw std::logic_error("Unimplemented DensityType in DensityTotal");
		}
		// This monitor pools between variables so ignores the dim it is passed:
		_dim.push_back(1);
    }

    void DensityTotal::update()
    {
		for (unsigned int ch = 0; ch < _nchain; ++ch) {
			double total = 0.0;
		    for (unsigned int i = 0; i < _nodes.size(); ++i) {
				total += _nodes[i]->logDensity(ch, PDF_FULL);
		    }
			if (total == JAGS_NA) {
			    // Don't try and convert NA to density or deviance
			}else if( _density_type == DENSITY ) {
				total = std::exp(total);
			}
			else if ( _density_type == DEVIANCE ) {
				total = -2.0 * total;
			}
		    _values[ch].push_back(total);
		}
    }
	
    vector<double> const &DensityTotal::value(unsigned int chain) const
    {
	return _values[chain];
    }

    vector<unsigned int> DensityTotal::dim() const
    {
	return _dim;
    }

    bool DensityTotal::poolChains() const
    {
	return false;
    }

    bool DensityTotal::poolIterations() const
    {
	return false;
    }

}}
