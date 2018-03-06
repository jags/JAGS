#include <config.h>

#include "PenaltyPV.h"

#include <util/nainf.h>
#include <distribution/Distribution.h>
// Required for PDFtype enum

using std::vector;
using std::string;

namespace jags {
namespace dic {

	// Public constructor:
    PenaltyPV::PenaltyPV(vector<Node const *> const &nodes, 
		vector<unsigned int> dim, string const &monitor_name)
	: Monitor(monitor_name, nodes), _nodes(nodes),
		_nchain(nodes[0]->nchain()), 
		_n(0), _pv(1,0.0), _mm(0.0), _mean(0.0), 
		_dim(vector<unsigned int> (1,1))
    {
		// This monitor pools between variables so ignores the dim it is passed
	}

    void PenaltyPV::update()
    {

	    for (unsigned int ch = 0; ch < _nchain; ++ch) {
			double newval = 0.0;
			for (unsigned int i = 0; i < _nodes.size(); ++i) {
				newval += (-2.0 * _nodes[i]->logDensity(ch, PDF_FULL));
			}
			if (newval == JAGS_NA) {
			    _mean = JAGS_NA;
			    _mm = JAGS_NA;
			    _pv[1] = JAGS_NA;
			}
			else {
				// _n is incremented after the loop so the current sample is:
				unsigned int effn = _n + ch + 1;
				double delta = newval - _mean;
				
				_mean += delta / effn;
				_mm += delta * (newval - _mean);
			}			
		}

		// Only necessary to update this once per iteration:
		_pv[0] = _mm / ( 2.0 * (double) (_n + _nchain - 1) );
		
		// Here _n is the total number of iterations * chains:
  	  	_n += _nchain;
		// Note: it is incremented after the loop as each chain counts as an iteration
		
   }
	
    vector<double> const &PenaltyPV::value(unsigned int chain) const
    {
	return _pv;
    }

    vector<unsigned int> PenaltyPV::dim() const
    {
	return _dim;
    }

    bool PenaltyPV::poolChains() const
    {
	return true;
    }

    bool PenaltyPV::poolIterations() const
    {
	return true;
    }

}}
