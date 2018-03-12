#include <config.h>
#include <graph/Node.h>
#include <distribution/Distribution.h>
// Required for PDFtype enum
#include <util/nainf.h>

#include "DensityMean.h"

#include <cmath>
#include <stdexcept>

using std::vector;
using std::string;
using std::exp;
using std::logic_error;

namespace jags {
namespace dic {

    DensityMean::DensityMean(vector<Node const *> const &nodes, vector<unsigned int> dim,
		DensityType const density_type, string const &monitor_name)
	: Monitor(monitor_name, nodes), _nodes(nodes),
	  _values(nodes[0]->nchain(), vector<double>(nodes.size(), 0.0)),
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
			throw logic_error("Unimplemented DensityType in DensityMean");
		}
		if ( monitor_name.compare(0, cdt.length(), cdt) != 0 ) {
			throw logic_error("Incorrect density type reported in monitor_name for DensityMean");
		}
		
		if ( monitor_name.find("_mean") == string::npos) {
			throw logic_error("Incorrect monitor type reported in monitor_name for DensityMean");
		}
    }

    void DensityMean::update()
    {
		_n++;
		for (unsigned int ch = 0; ch < _nchain; ++ch) {
		    vector<double> &rmean  = _values[ch];			
		    for (unsigned int i = 0; i < _nodes.size(); ++i) {
				double newval = _nodes[i]->logDensity(ch, PDF_FULL);
				if (newval == JAGS_NA) {
				    rmean[i] = JAGS_NA;
				}
				else {
					if( _density_type == DENSITY ) {
						newval = exp(newval);
					}
					else if ( _density_type == DEVIANCE ) {
						newval = -2.0 * newval;
					}
				    rmean[i] -= (rmean[i] - newval)/_n;
				}
		    }
		}
    }
	
    vector<double> const &DensityMean::value(unsigned int chain) const
    {
	return _values[chain];
    }

    vector<unsigned int> DensityMean::dim() const
    {
	return _dim;
    }

    bool DensityMean::poolChains() const
    {
	return false;
    }

    bool DensityMean::poolIterations() const
    {
	return true;
    }

}}
