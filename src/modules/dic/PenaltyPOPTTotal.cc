#include <config.h>

#include "PenaltyPOPTTotal.h"
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>

#include <cmath>

using std::vector;
using std::string;

namespace jags {

namespace dic {

	// Public constructor:
    PenaltyPOPTTotal::PenaltyPOPTTotal(vector<Node const *> const &nodes,
			 vector<unsigned int> dim, 
			 string const &monitor_name,
			 vector<RNG *> const &rngs,
			 unsigned int nrep)
	: PenaltyPDTotal(nodes, dim, monitor_name, rngs, nrep, 2.0),
			 _weights(nodes.size(), 0.0), _n(0)
    {
    }
	
    void PenaltyPOPTTotal::update()
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
			// Here this is the average weight (is sum for PenaltyPOPT):
			_weights[k] -= (_weights[k] - wsum) / _n;
			popt += pdsum / _weights[k];
		}
		
		popt *= _scale_cst;
		_values.push_back(popt);				
    }
	
	PenaltyPOPTTotal::~PenaltyPOPTTotal()
	{
	}
	

}}
