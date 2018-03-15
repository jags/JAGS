#include <config.h>
#include <graph/Node.h>
#include <distribution/Distribution.h>
// Required for PDFtype enum
#include <util/nainf.h>

#include "DensityTrace.h"

#include <cmath>
#include <stdexcept>

using std::vector;
using std::string;
using std::exp;
using std::logic_error;

namespace jags {
namespace dic {

    DensityTrace::DensityTrace(vector<Node const *> const &nodes, vector<unsigned int> const &dim,
		DensityType const density_type, string const &monitor_name)
	: Monitor(monitor_name, nodes), _nodes(nodes), _values(nodes[0]->nchain()),
	  _density_type(density_type), _dim(dim), _nchain(nodes[0]->nchain())

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
			throw logic_error("Unimplemented DensityType in DensityTrace");
		}
		if ( monitor_name.compare(0, cdt.length(), cdt) != 0 ) {
			throw logic_error("Incorrect density type reported in monitor_name for DensityTrace");
		}
		
		if ( monitor_name.find("_trace") == string::npos) {
			throw logic_error("Incorrect monitor type reported in monitor_name for DensityTrace");
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
					newval = exp(newval);
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
