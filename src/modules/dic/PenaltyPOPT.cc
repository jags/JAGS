#include <config.h>

#include "PenaltyPOPT.h"
#include <module/ModuleError.h>

#include <cmath>

using std::vector;
using std::string;

namespace jags {
namespace dic {


	PenaltyPOPT::PenaltyPOPT(vector<Node const *> const &nodes,
			 vector<unsigned int> dim, 
			 string const &monitor_name,
			 vector<RNG *> const &rngs,
			 unsigned int nrep)
		: PenaltyPD(nodes, dim, monitor_name, rngs, nrep, 2.0),
			 _weights(nodes.size(), 0)
    {
		if (_nchain < 2) {
		    throwLogicError("The popt monitor needs at least 2 chains");
		}
    }
	
    void PenaltyPOPT::update()
    {
		// Not actually needed for popt (just for pD):
		// _n++;
		
		vector<double> w(_nchain);
		for (unsigned int k = 0; k < _values.size(); ++k) {
	    
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
		
		    pdsum /= wsum;
		    pdsum *= _scale_cst;

		    _weights[k] += wsum;
		    _values[k] += wsum * (pdsum - _values[k])/_weights[k];
			
		}
	}

}}
