#include <module/Module.h>

#include "DevianceMonitorFactory.h"
#include "PDMonitorFactory.h"
#include "PDTraceFactory.h"
//#include "WAICMonitorFactory.h"
#include "NodeDensityMonitorFactory.h"
#include "ObsStochDensMonitorFactory.h"

using std::vector;

namespace jags {
namespace dic {

    class DICModule: public Module {
    public:
	DICModule();
	~DICModule();
    };
    
    DICModule::DICModule() 
	: Module("dic")
    {
	
		/*  DevianceMonitorFactory (and DevianceTrace/DevianceMean) from 
			JAGS 4.3.0 could now be retired?
			The same quantities are provided by ObsStochDensMonitorFactory 
			and are named for backwards-compatibility (althouh currently
			disabled to avoid clashes)
		*/
		insert(new DevianceMonitorFactory);

		// Handles density-related monitors for a given node:
		insert(new NodeDensityMonitorFactory);
		// Handles density-related monitors for all observed stochastic nodes:
		insert(new ObsStochDensMonitorFactory);

		// Unchanged from JAGS 4.3.0
		insert(new PDMonitorFactory);
		// Unchanged from JAGS 4.3.0
		insert(new PDTraceFactory);

		// Now retired:
	//	insert(new WAICMonitorFactory);
    }
    
    DICModule::~DICModule() {
	
	vector<MonitorFactory*> const &mvec = monitorFactories();
	for (unsigned int i = 0; i < mvec.size(); ++i) {
	    delete mvec[i];
	}
    }

}}

jags::dic::DICModule _dic_module;

