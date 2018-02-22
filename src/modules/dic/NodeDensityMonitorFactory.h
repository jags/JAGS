#ifndef NODE_DENSITY_MONITOR_FACTORY_H_
#define NODE_DENSITY_MONITOR_FACTORY_H_

#include <model/MonitorFactory.h>
#include <model/Monitor.h>
// Required for DensityType enum:
#include <model/NodeArraySubset.h>

namespace jags {
namespace dic {

	enum MonitorType {TRACE, MEAN, VARIANCE, TOTAL};

    class NodeDensityMonitorFactory : public MonitorFactory
    {
      public:
	Monitor *getMonitor(std::string const &name, Range const &range,
			    BUGSModel *model, std::string const &type,
			    std::string &msg);
	std::string name() const;
    };
    
}}

#endif /* NODE_DENSITY_MONITOR_FACTORY_H_ */
