#ifndef LOGDENSITY_MONITOR_FACTORY_H_
#define LOGDENSITY_MONITOR_FACTORY_H_

#include <model/MonitorFactory.h>

namespace jags {
namespace base {

    class LogDensityMonitorFactory : public MonitorFactory
    {
      public:
	Monitor *getMonitor(std::string const &name, Range const &range, 
			    BUGSModel *model, std::string const &type,
			    std::string &msg);
	std::string name() const;
    };
    
}}

#endif /* LOGDENSITY_MONITOR_FACTORY_H_ */
