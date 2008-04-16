#ifndef POPT_MONITOR_H_
#define POPT_MONITOR_H_

#include <model/Monitor.h>
#include <graph/StochasticNode.h>

#include <vector>

class StochasticNode;
class RNG;

namespace dic {

    class PoptMonitor : public Monitor {
    protected:
        StochasticNode const *_snode;
	std::vector<double> _weights; // weights for each chain
	std::vector<double> _values; // sampled values
    public:
	PoptMonitor(StochasticNode const *snode,
		    unsigned int start,  unsigned int thin);
	unsigned int nchain() const;
	std::vector<unsigned int> dim() const;
	std::vector<double> const &value(unsigned int chain) const;
	void reserve(unsigned int niter);
	SArray dump() const;
    };

}

#endif /* POPT_MONITOR_H_ */
