#ifndef DEFAULT_PD_MONITOR_H_
#define DEFAULT_PD_MONITOR_H_

#include "PDMonitor.h"

#include <model/Monitor.h>
#include <graph/StochasticNode.h>
#include <vector>

class StochasticNode;
class RNG;

namespace dic {

    class DefaultPDMonitor : public PDMonitor {
	StochasticNode _repnode;
	const std::vector<RNG *> _rngs;
	unsigned int _nrep;
    public:
	DefaultPDMonitor(StochasticNode const *snode,
			 unsigned int start,  unsigned int thin, 
			 std::vector<RNG *> const &rngs, unsigned int nrep); 
	void doUpdate();
    };

}

#endif /* DEFAULT_PD_MONITOR_H_ */
