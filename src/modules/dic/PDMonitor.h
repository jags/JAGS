#ifndef PD_MONITOR_H_
#define PD_MONITOR_H_

#include <model/Monitor.h>
#include <graph/StochasticNode.h>

#include <vector>

class StochasticNode;
class RNG;

namespace dic {

    class PDMonitor : public Monitor {
	StochasticNode _repnode;
	const std::vector<RNG *> _rngs;
	unsigned int _nrep;
	std::vector<double> _values; // sampled values
    public:
	PDMonitor(StochasticNode const *snode,
		  unsigned int start,  unsigned int thin, 
		  std::vector<RNG *> const &rngs, unsigned int nrep); 
	unsigned int nchain() const;
	std::vector<unsigned int> dim() const;
	std::vector<double> const &value(unsigned int chain) const;
	void doUpdate();
	void reserve(unsigned int niter);
	SArray dump() const;
    };

}

#endif /* PD_MONITOR_H_ */
