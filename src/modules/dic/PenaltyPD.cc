#include <config.h>

#include "PenaltyPD.h"
#include <module/ModuleError.h>

using std::vector;
using std::string;

namespace jags {

namespace dic {

	// Public constructor:
    PenaltyPD::PenaltyPD(vector<Node const *> const &nodes,
			 vector<unsigned long> const &dim, 
			 string const &monitor_name,
			 vector<RNG *> const &rngs, unsigned int nrep)
	: Monitor(monitor_name, nodes), _nodes(nodes), _rngs(rngs),
	  _nrep(nrep), _values(nodes.size(), 0.0), _dim(dim), _scale_cst(0.5),
	  _nchain(rngs.size()), _n(0)

    {
		if (_nchain < 2) {
		    throwLogicError("The pD monitor needs at least 2 chains");
		}
    }

	// Protected constructor which PenaltyPOPT uses:
    PenaltyPD::PenaltyPD(vector<Node const *> const &nodes,
			 vector<unsigned long> const &dim, 
			 string const &monitor_name,
			 vector<RNG *> const &rngs,
			 unsigned int nrep, double scale)
	: Monitor(monitor_name, nodes), _nodes(nodes), _rngs(rngs),
	  _nrep(nrep),_values(nodes.size(), 0.0), _dim(dim), _scale_cst(1.0),
	  _nchain(rngs.size()), _n(0)

    {
		if (_nchain < 2) {
		    throwLogicError("The popt monitor needs at least 2 chains");
		}
    }


    PenaltyPD::~PenaltyPD() 
    {
    }

    vector<unsigned long> PenaltyPD::dim() const
    {
	return _dim;
    }
 
    vector<double> const &PenaltyPD::value(unsigned int chain) const
    {
	return _values;
    }

    bool PenaltyPD::poolChains() const
    {
	return true;
    }

    bool PenaltyPD::poolIterations() const
    {
	return true;
    }

    void PenaltyPD::update()
    {
		_n++;
		for (unsigned int k = 0; k < _values.size(); ++k) {
	    
		    double pdsum = 0;
		    for (unsigned int i = 0; i < _nchain; ++i) {
				for (unsigned int j = 0; j < i; ++j) {
				    pdsum += _nodes[k]->KL(i, j, _rngs[i], _nrep);
					pdsum += _nodes[k]->KL(j, i, _rngs[j], _nrep);
				}
		    }
			// Number of combinations of chains * unity weight products:
			// pdsum /= (double) _nchain * (_nchain - 1) * 0.5;
		    // pdsum *= _scale_cst;
			// i.e. equivalent to:
			pdsum /= (double) _nchain * (_nchain - 1);

		    _values[k] -= (_values[k] - pdsum)/_n;
		
		}
    }

}}
