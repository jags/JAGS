#ifndef PENALTY_POPT_TOTAL_H
#define PENALTY_POPT_TOTAL_H

#include <model/Monitor.h>
#include <graph/Node.h>
#include <rng/RNG.h>

#include "PenaltyPDTotal.h"

#include <vector>

/* 	This version of the POPTTotal monitor is far more efficient
	but the first few estimates of POPT will use a poorly estimated
	average weight.  Also the trace mean won't exactly correspond to the
	sum of the POPT monitor (although it is close except over a short
	run)
*/

namespace jags {
namespace dic {

    class PenaltyPOPTTotal : public PenaltyPDTotal {
		unsigned int _n;
		std::vector<double> _weights;
    public:
	PenaltyPOPTTotal(std::vector<Node const *> const &nodes,
		  std::vector<unsigned int> dim, 
		  std::string const &monitor_name,
		  std::vector<RNG *> const &rngs,
		  unsigned int nrep);

	void update();
	~PenaltyPOPTTotal();
    };

}}

#endif /* PENALTY_POPT_TOTAL_H */
