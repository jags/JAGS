#ifndef WAIC_MONITOR_FACTORY_H_
#define WAIC_MONITOR_FACTORY_H_

#include <model/MonitorFactory.h>

namespace jags {
    namespace dic {

	class WAICMonitorFactory : public MonitorFactory
	{
	public:
	    Monitor *getMonitor(std::string const &name, Range const &range,
				BUGSModel *model, std::string const &type,
				std::string &msg);
	    std::string name() const;
	};
    
    }
}

#endif /* WAIC_MONITOR_FACTORY_H_ */
