#include <config.h>

#include "PenaltyPOPTTotalRep.h"
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>

#include <cmath>

using std::vector;
using std::string;

namespace jags {
namespace dic {

	// Public constructor:
    PenaltyPOPTTotalRep::PenaltyPOPTTotalRep(vector<Node const *> const &nodes,
			 string const &monitor_name,
			 vector<RNG *> const &rngs,
			 unsigned int nrep)
	: PenaltyPDTotal(nodes, monitor_name, rngs, nrep, 2.0),
	  _n(0), _weights(nodes.size(), 0.0), _nodetrace(nodes.size())
    {
		/* This is a hack to allow the total popt to be adjusted by
		the running mean weight when requested by const value() */
		_totalpopt = new vector<double>;
    }
	
    void PenaltyPOPTTotalRep::update()
    {
		_n++;

		double popt = 0.0;
		
		vector<double> w(_nchain);
		for (unsigned int k = 0; k < _nodes.size(); ++k) {
			
		    double pdsum = 0;
		    double wsum = 0;
		    for (unsigned int i = 0; i < _nchain; ++i) {
				w[i] = std::exp(- _nodes[k]->logDensity(i, PDF_FULL));
				for (unsigned int j = 0; j < i; ++j) {
				    pdsum += w[i] * w[j] * (
					_nodes[k]->KL(i, j, _rngs[i], _nrep) +
					_nodes[k]->KL(j, i, _rngs[j], _nrep));
				    wsum += w[i] * w[j];
				}
		    }
			
			_nodetrace[k].push_back(pdsum);
			
			// Here this is the average weight (is sum for PenaltyPOPT):
			_weights[k] -= (_weights[k] - wsum) / _n;
			popt += pdsum / _weights[k];
		}
    }
	
    vector<double> const &PenaltyPOPTTotalRep::value(unsigned int chain) const
    {
		// Adjust by the running mean weights:
		(*_totalpopt).resize(_n);
		for(unsigned int i = 0; i < _n; ++i){
			(*_totalpopt)[i] = 0.0;
			for(unsigned int k = 0; k < _nodes.size(); ++k){
				(*_totalpopt)[i] += (_nodetrace[k][i] / _weights[k]);
			}
			(*_totalpopt)[i] *= _scale_cst;
		}
		
		return (*_totalpopt);
		
		// Note: _values isn't actually used anywhere
		// return _values;
    }

	PenaltyPOPTTotalRep::~PenaltyPOPTTotalRep()
	{
		// Prevent memory leak with the hack:
		delete _totalpopt;
		_totalpopt = 0;
	}
	

}}
