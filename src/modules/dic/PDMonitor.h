#ifndef PD_MONITOR_H_
#define PD_MONITOR_H_

#include <model/Monitor.h>
#include <graph/StochasticNode.h>

#include <vector>

class StochasticNode;
class RNG;

namespace dic {

    class PDMonitor : public Monitor {
    protected:
	std::vector<double> _values;
    public:
	PDMonitor(StochasticNode const *node, unsigned int start,
		  unsigned int thin);
	unsigned int nchain() const;
	std::vector<unsigned int> dim() const;
	std::vector<double> const &value(unsigned int chain) const;
	void reserve(unsigned int niter);
	SArray dump() const;
    };

}

#endif /* PD_MONITOR_H_ */
