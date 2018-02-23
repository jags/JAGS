#include <config.h>
#include <graph/Node.h>
#include <distribution/Distribution.h>
// Required for PDFtype enum
#include <util/nainf.h>

#include "DensityVariance.h"

#include <cmath>

using std::vector;
using std::string;

namespace jags {
namespace dic {

    DensityVariance::DensityVariance(vector<Node const *> const &nodes, vector<unsigned int> dim,
		DensityType const density_type, string const &monitor_name)
	: Monitor(monitor_name, nodes), _nodes(nodes), _density_type(density_type), 
		_nchain(nodes[0]->nchain()), 
		_means(nodes[0]->nchain(), vector<double>(nodes.size(), 0.0)),
		_mms(nodes[0]->nchain(), vector<double>(nodes.size(), 0.0)),
		_variances(nodes[0]->nchain(), vector<double>(nodes.size(), 0.0)), _n(0), _dim(dim)
    {
		if( _density_type != DENSITY && _density_type != LOGDENSITY && _density_type != DEVIANCE ) {
			throw std::logic_error("Unimplemented DensityType in DensityVariance");
		}
    }

    void DensityVariance::update()
    {
		_n++;
		for (unsigned int ch = 0; ch < _nchain; ++ch) {
		    vector<double> &rmean  = _means[ch];
		    vector<double> &rmm  = _mms[ch];
			vector<double> &rvar  = _variances[ch];		
		    for (unsigned int i = 0; i < _nodes.size(); ++i) {
				double newval = _nodes[i]->logDensity(ch, PDF_FULL);
				if (newval == JAGS_NA) {
				    rmean[i] = JAGS_NA;
					rmm[i] = JAGS_NA;
					rvar[i] = JAGS_NA;
				}
				else {
					if( _density_type == DENSITY ) {
						newval = std::exp(newval);
					}
					else if ( _density_type == DEVIANCE ) {
						newval = -2.0 * newval;
					}

					double delta = newval - rmean[i];
					rmean[i] += delta / _n;
					rmm[i] += delta * (newval - rmean[i]);
					rvar[i] = rmm[i] / (double) (_n - 1);
				}
		    }
		}
    }
	
    vector<double> const &DensityVariance::value(unsigned int chain) const
    {
	return _variances[chain];
    }

    vector<unsigned int> DensityVariance::dim() const
    {
	return _dim;
    }

    bool DensityVariance::poolChains() const
    {
	return false;
    }

    bool DensityVariance::poolIterations() const
    {
	return true;
    }

}}
