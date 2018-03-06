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
		
		// Last monitor factories to be used are given first:
		// density-related monitors for a given node:
		insert(new NodeDensityMonitorFactory);
		// density-related monitors for all observed stochastic nodes:
		insert(new ObsStochDensMonitorFactory);

		/*  The below factories from JAGS 4.3.0 could now be retired
			The same quantities are provided by ObsStochDensMonitorFactory 
			and are named for backwards-compatibility (although currently
			disabled to avoid clashes)
		*/
		
		// To compile the new equivalents use the flag -D NEW_DIC
		#ifndef NEW_DIC
			insert(new DevianceMonitorFactory);

			// Unchanged from JAGS 4.3.0
			insert(new PDMonitorFactory);
			// Unchanged from JAGS 4.3.0
			insert(new PDTraceFactory);
		#else
			printf("  [Using new penalty monitors]  \n");
							
		#endif
		// Otherwise the old versions remain and take precedence
		// TODO: remove this macro and printf
				

		// Now retired:
		// insert(new WAICMonitorFactory);
	
    }
    
    DICModule::~DICModule() {
	
	vector<MonitorFactory*> const &mvec = monitorFactories();
	for (unsigned int i = 0; i < mvec.size(); ++i) {
	    delete mvec[i];
	}
    }

}}

jags::dic::DICModule _dic_module;

