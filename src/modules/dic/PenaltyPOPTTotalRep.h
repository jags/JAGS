#ifndef PENALTY_POPT_TOTAL_REP_H
#define PENALTY_POPT_TOTAL_REP_H

#include <model/Monitor.h>
#include <graph/Node.h>
#include <rng/RNG.h>

#include "PenaltyPDTotal.h"

#include <vector>

/* 	This version of the POPTTotal monitor is far less efficient
	and also includes a bit of a hack to allow the const value
	method to trigger calculations based on the long-run average
	weights.  But the trace mean will exactly correspond to the
	sum of the POPT monitor.
*/

namespace jags {
namespace dic {

    class PenaltyPOPTTotalRep : public PenaltyPDTotal {
		unsigned int _n;
		std::vector<double> _weights;
		/* This is a hack to allow the total popt to be adjusted by
		the running mean weight when requested by const value() */
		std::vector<double>* _totalpopt;
		std::vector< std::vector<double> > _nodetrace;
    public:
	PenaltyPOPTTotalRep(std::vector<Node const *> const &nodes,
		  std::vector<unsigned int> dim, 
		  std::string const &monitor_name,
		  std::vector<RNG *> const &rngs,
		  unsigned int nrep);

  	std::vector<double> const &value(unsigned int chain) const;
	void update();
	~PenaltyPOPTTotalRep();
    };

}}

#endif /* PENALTY_POPT_TOTAL_REP_H */
