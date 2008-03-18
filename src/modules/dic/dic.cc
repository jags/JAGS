#include <Module.h>

#include "DevianceMonitorFactory.h"
#include "PDMonitorFactory.h"
#include "PoptMonitorFactory.h"

using std::vector;

namespace dic {

class DICModule: public Module {
public:
    DICModule();
    ~DICModule();
};

DICModule::DICModule() {

    insert(new DevianceMonitorFactory);
    insert(new PDMonitorFactory);
    insert(new PoptMonitorFactory);
}

DICModule::~DICModule() {

    vector<MonitorFactory*> const &mvec = monitorFactories();
    for (unsigned int i = 0; i < mvec.size(); ++i) {
	delete mvec[i];
    }
}

}

dic::DICModule _dic_module;

