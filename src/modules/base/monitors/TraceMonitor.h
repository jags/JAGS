#ifndef TRACE_MONITOR_H_
#define TRACE_MONITOR_H_

#include <model/Monitor.h>
#include <model/NodeArraySubset.h>

#include <vector>

namespace jags {
    namespace base {

	/**
	 * @short Stores sampled values of a given Node
	 */
	class TraceMonitor : public Monitor {
	    NodeArraySubset _subset;
	    std::vector<std::vector<double> > _values; // sampled values
	  public:
	    TraceMonitor(NodeArraySubset const &subset);
	    void update();
	    std::vector<double> const &value(unsigned int chain) const;
	    std::vector<unsigned long> dim() const;
	    bool poolChains() const;
	    bool poolIterations() const;
	};
	
    }
}

#endif /* TRACE_MONITOR_H_ */
