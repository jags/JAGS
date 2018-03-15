#include <config.h>
#include <graph/Node.h>
#include <distribution/Distribution.h>
// Required for PDFtype enum
#include <util/nainf.h>

#include "DensityTotal.h"

#include <cmath>
#include <stdexcept>

using std::vector;
using std::string;
using std::logic_error;
using std::exp;

namespace jags {
namespace dic {

    DensityTotal::DensityTotal(vector<Node const *> const &nodes, vector<unsigned int> const &dim,
		DensityType const density_type, string const &monitor_name)
	: Monitor(monitor_name, nodes), _nodes(nodes), _values(nodes[0]->nchain()),
	  _density_type(density_type), _dim(vector<unsigned int> (1,1)), _nchain(nodes[0]->nchain())
    {
		// This monitor pools between variables so ignores the dim it is passed

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
			throw logic_error("Unimplemented DensityType in DensityTotal");
		}

		// Required for back-compatibility (only from ObsStochDensMonitorFactory):
		if ( monitor_name == "trace" ) {
			if ( _density_type != DEVIANCE ) {
				throw logic_error("DensityTotal is reporting a non-DEVIANCE type with monitor_name trace");
			}
		}
		else {

			if ( monitor_name.compare(0, cdt.length(), cdt) != 0 ) {
				throw logic_error("Incorrect density type reported in monitor_name for DensityTotal");
			}
		
			if ( monitor_name.find("_total") == string::npos) {
				throw logic_error("Incorrect monitor type reported in monitor_name for DensityTotal");
			}
		
		}

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
				total = exp(total);
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
