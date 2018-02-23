#include <config.h>
#include <graph/Node.h>
#include <distribution/Distribution.h>
// Required for PDFtype enum
#include <util/nainf.h>

#include "DensityPoolMean.h"

#include <cmath>

using std::vector;
using std::string;

namespace jags {
namespace dic {

    DensityPoolMean::DensityPoolMean(vector<Node const *> const &nodes,
		DensityType const density_type, string const &monitor_name)
	: Monitor(monitor_name, nodes), _nodes(nodes), _density_type(density_type), 
		_nchain(nodes[0]->nchain()), _values(nodes.size(), 0.0), _n(0)
    {
		if( _density_type != DENSITY && _density_type != LOGDENSITY && _density_type != DEVIANCE ) {
			throw std::logic_error("Unimplemented DensityType in DensityPoolMean");
		}
    }

    void DensityPoolMean::update()
    {
		_n++;
		for (unsigned int i = 0; i < _nodes.size(); ++i) {
		    double newval = 0.0;
		    for (unsigned int ch = 0; ch < _nchain; ++ch) {
				newval += _nodes[i]->logDensity(ch, PDF_FULL) / _nchain;
		    }
			if (newval == JAGS_NA) {
			    _values[i] = JAGS_NA;
			}
			else {
				if( _density_type == DENSITY ) {
					newval = std::exp(newval);
				}
				else if ( _density_type == DEVIANCE ) {
					newval = -2.0 * newval;
				}
			    _values[i] -= (_values[i] - newval)/_n;
			}
		}
    }
	
    vector<double> const &DensityPoolMean::value(unsigned int chain) const
    {
	return _values;
    }

    vector<unsigned int> DensityPoolMean::dim() const
    {
	return vector<unsigned int>(1, _nodes.size());
    }

    bool DensityPoolMean::poolChains() const
    {
	return true;
    }

    bool DensityPoolMean::poolIterations() const
    {
	return true;
    }

}}
