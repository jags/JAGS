#include <config.h>
#include <graph/Node.h>
#include <distribution/Distribution.h>
// Required for PDFtype enum
#include <util/nainf.h>

#include "DensityPoolMean.h"

#include <cmath>
#include <stdexcept>

using std::vector;
using std::string;
using std::exp;
using std::logic_error;

namespace jags {
namespace dic {

    DensityPoolMean::DensityPoolMean(vector<Node const *> const &nodes, vector<unsigned int> dim,
		DensityType const density_type, string const &monitor_name)
	: Monitor(monitor_name, nodes), _nodes(nodes), _values(nodes.size(), 0.0),
	  _density_type(density_type), _dim(dim), _nchain(nodes[0]->nchain()), _n(0) 
    {
		// Sanity check that input arguments match to this function:
		
		string cdt("nomatch");
		if ( _density_type == DENSITY ) {
			cdt.assign("density");			
		}
		else if ( _density_type == LOGDENSITY ) {
			cdt.assign("logdensity");			
		}
		else if ( _density_type == DEVIANCE ) {
			cdt.assign("deviance");			
		}
		else {
			throw logic_error("Unimplemented DensityType in DensityPoolMean");
		}
		
		// Required for back-compatibility (only from ObsStochDensMonitorFactory):
		if ( monitor_name == "mean" ) {
			if ( _density_type != DEVIANCE ) {
				throw logic_error("DensityPoolMean is reporting a non-DEVIANCE type with monitor_name mean");
			}
		}
		else {
				
			if ( monitor_name.compare(0, cdt.length(), cdt) != 0 ) {
				throw logic_error("Incorrect density type reported in monitor_name for DensityPoolMean");
			}		
			if ( monitor_name.find("_poolmean") == string::npos) {
				throw logic_error("Incorrect monitor type reported in monitor_name for DensityPoolMean");
			}
			
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
					newval = exp(newval);
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
	return _dim;
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
