#ifndef DEFAULT_POPT_MONITOR_H_
#define DEFAULT_POPT_MONITOR_H_

#include <model/Monitor.h>
#include <graph/StochasticNode.h>

#include <vector>

class StochasticNode;
class RNG;

namespace dic {

    class PoptMonitor : public Monitor {
	StochasticNode const *_snode;
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
