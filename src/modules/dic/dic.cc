#include <module/Module.h>

#include "DevianceMonitorFactory.h"
#include "PDMonitorFactory.h"
#include "PDTraceFactory.h"
#include "WAICMonitorFactory.h"
#include "NodeDensityMonitorFactory.h"

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
	
	insert(new DevianceMonitorFactory);
	insert(new PDMonitorFactory);
	insert(new PDTraceFactory);
	insert(new WAICMonitorFactory);
	insert(new NodeDensityMonitorFactory);
    }
    
    DICModule::~DICModule() {
	
	vector<MonitorFactory*> const &mvec = monitorFactories();
	for (unsigned int i = 0; i < mvec.size(); ++i) {
	    delete mvec[i];
	}
    }

}}

jags::dic::DICModule _dic_module;

