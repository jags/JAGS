#ifndef LOGDENSITY_MONITOR_H_
#define LOGDENSITY_MONITOR_H_

#include <model/Monitor.h>
#include <model/NodeArraySubset.h>

#include <vector>

namespace jags {
    namespace base {

	/**
	 * @short Stores log densities corresponding to sampled values of a given Node
	 */
	class LogDensityMonitor : public Monitor {
	    NodeArraySubset _subset;
	    std::vector<std::vector<double> > _values; // log densities corresponding to sampled values
	  public:
	    LogDensityMonitor(NodeArraySubset const &subset);
	    void update();
	    std::vector<double> const &value(unsigned int chain) const;
	    std::vector<unsigned int> dim() const;
	    bool poolChains() const;
	    bool poolIterations() const;
	};
	
    }
}

#endif /* LOGDENSITY_MONITOR_H_ */
