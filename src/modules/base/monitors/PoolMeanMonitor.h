#ifndef POOLMEAN_MONITOR_H_
#define POOLMEAN_MONITOR_H_

#include <model/Monitor.h>
#include <model/NodeArraySubset.h>

#include <vector>

namespace jags {
namespace base {

    /**
     * @short Stores running mean (pooled between chains) of a given Node
     */
    class PoolMeanMonitor : public Monitor {
	NodeArraySubset _subset;
	std::vector<double> _values; // sampled values
	unsigned int _n;
    public:
	PoolMeanMonitor(NodeArraySubset const &subset);
	void update();
	std::vector<double> const &value(unsigned int chain) const;
	std::vector<unsigned long> dim() const;
	bool poolChains() const;
	bool poolIterations() const;
    };

}}

#endif /* POOLMEAN_MONITOR_H_ */
