#ifndef PENALTY_POPT_H_
#define PENALTY_POPT_H_

#include "PenaltyPD.h"

#include <vector>

namespace jags {
namespace dic {

   class PenaltyPOPT : public PenaltyPD {
		std::vector<double> _weights;
    public:
	PenaltyPOPT(std::vector<Node const *> const &nodes,
		  std::vector<unsigned long> const &dim, 
		  std::string const &monitor_name,
		  std::vector<RNG *> const &rngs,
		  unsigned int nrep);

	void update();
    };

}}

#endif /* PENALTY_POPT_H_ */
