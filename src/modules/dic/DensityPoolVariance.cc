#include <config.h>
#include <graph/Node.h>
#include <distribution/Distribution.h>
// Required for PDFtype enum
#include <util/nainf.h>

#include "DensityPoolVariance.h"

#include <cmath>

using std::vector;
using std::string;

namespace jags {
namespace dic {

    DensityPoolVariance::DensityPoolVariance(vector<Node const *> const &nodes, vector<unsigned int> dim,
		DensityType const density_type, string const &monitor_name)
	: Monitor(monitor_name, nodes), _nodes(nodes), _density_type(density_type), 
		_nchain(nodes[0]->nchain()), 
		_means(nodes.size(), 0.0),
		_mms(nodes.size(), 0.0),
		_variances(nodes.size(), 0.0), _n(0), _dim(dim)
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
			throw std::logic_error("Unimplemented DensityType in DensityPoolVariance");
		}
		if ( monitor_name.compare(0, cdt.length(), cdt) != 0 ) {
			throw std::logic_error("Incorrect density type reported in monitor_name for DensityPoolVariance");
		}
		
		if ( monitor_name.find("_poolvariance") == string::npos) {
			throw std::logic_error("Incorrect monitor type reported in monitor_name for DensityPoolVariance");
		}
    }

    void DensityPoolVariance::update()
    {

		for (unsigned int i = 0; i < _nodes.size(); ++i) {
		    for (unsigned int ch = 0; ch < _nchain; ++ch) {
				double newval = _nodes[i]->logDensity(ch, PDF_FULL);
				if (newval == JAGS_NA) {
				    _means[i] = JAGS_NA;
				    _mms[i] = JAGS_NA;
				    _variances[i] = JAGS_NA;
				}
				else {
					if( _density_type == DENSITY ) {
						newval = std::exp(newval);
					}
					else if ( _density_type == DEVIANCE ) {
						newval = -2.0 * newval;
					}
					
					// _n is incremented after the loop so the current sample is:
					unsigned int effn = _n + ch + 1;
					double delta = newval - _means[i];
					
					_means[i] += delta / effn;
					_mms[i] += delta * (newval - _means[i]);
				}
			}
			
			// Only necessary to update this once per iteration:
			_variances[i] = _mms[i] / (double) (_n + _nchain - 1);			
		}
		
		// Here _n is the total number of iterations * chains:
  	  	_n += _nchain;
		// Note: it is incremented after the loop as each chain counts as an iteration
		
   }
	
    vector<double> const &DensityPoolVariance::value(unsigned int chain) const
    {
	return _variances;
    }

    vector<unsigned int> DensityPoolVariance::dim() const
    {
	return _dim;
    }

    bool DensityPoolVariance::poolChains() const
    {
	return true;
    }

    bool DensityPoolVariance::poolIterations() const
    {
	return true;
    }

}}
