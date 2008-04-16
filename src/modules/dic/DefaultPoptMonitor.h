#ifndef DEFAULT_POPT_MONITOR_H_
#define DEFAULT_POPT_MONITOR_H_

#include "PoptMonitor.h"

class StochasticNode;
class RNG;

namespace dic {

    class DefaultPoptMonitor : public PoptMonitor {
	StochasticNode _repnode;
	const std::vector<RNG *> _rngs;
	unsigned int _nrep;
    public:
	DefaultPoptMonitor(StochasticNode const *snode,
			   unsigned int start,  unsigned int thin, 
			   std::vector<RNG *> const &rngs, unsigned int nrep); 
	void doUpdate();
    };

}

#endif /* DEFAULT_POPT_MONITOR_H_ */
