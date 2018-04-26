#include <config.h>

#include "PenaltyPDTotal.h"
#include <module/ModuleError.h>

using std::vector;
using std::string;

namespace jags {
namespace dic {

	// Public constructor:
    PenaltyPDTotal::PenaltyPDTotal(vector<Node const *> const &nodes,
			 string const &monitor_name,
			 vector<RNG *> const &rngs, unsigned int nrep)
	: Monitor(monitor_name, nodes), _nodes(nodes), _rngs(rngs),
	  _nrep(nrep), _nchain(rngs.size()), _values(),
	  _dim(vector<unsigned long> (1,1)), 
	  _scale_cst(1.0/2.0)
    {
		// This monitor pools between variables so ignores the dim it is passed
		if (_nchain < 2) {
		    throwLogicError("The pD total monitor needs at least 2 chains");
		}
    }

	// Protected constructor which PenaltyPOPTTotal uses:
    PenaltyPDTotal::PenaltyPDTotal(vector<Node const *> const &nodes,
			 string const &monitor_name,
			 vector<RNG *> const &rngs,
			 unsigned int nrep, double scale)
	: Monitor(monitor_name, nodes), _nodes(nodes), _rngs(rngs),
	  _nrep(nrep), _nchain(rngs.size()), _values(),
	  _dim(vector<unsigned long> (1,1)), _scale_cst(scale/2.0)
    {
		// This monitor pools between variables so ignores the dim it is passed
		if (_nchain < 2) {
		    throwLogicError("The popt total monitor needs at least 2 chains");
		}
    }

    PenaltyPDTotal::~PenaltyPDTotal() 
    {
    }

    vector<unsigned long> PenaltyPDTotal::dim() const
    {
	return _dim;
    }
 
    vector<double> const &PenaltyPDTotal::value(unsigned int chain) const
    {
	return _values;
    }

    bool PenaltyPDTotal::poolChains() const
    {
	return true;
    }

    bool PenaltyPDTotal::poolIterations() const
    {
	return false;
    }

    void PenaltyPDTotal::update()
    {

	double pd = 0;
	for (unsigned int k = 0; k < _nodes.size(); ++k) {
	    for (unsigned int i = 0; i < _nchain; ++i) {
		for (unsigned int j = 0; j < i; ++j) {
		    pd += _nodes[k]->KL(i, j, _rngs[i], _nrep);
		    pd += _nodes[k]->KL(j, i, _rngs[j], _nrep);
		}
	    }
	}

	// NB: constant multiplier 2/_scale_cst removed:
	pd /= _nchain * (_nchain - 1);
			
	_values.push_back(pd);
    }

}}
