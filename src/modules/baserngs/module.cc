#include <Module.h>

#include "BaseRNGFactory.h"

using std::vector;

class BaseRNGModule : public Module {
public:
    BaseRNGModule();
    ~BaseRNGModule();
};

BaseRNGModule::BaseRNGModule() {

    insert(new BaseRNGFactory);
}

BaseRNGModule::~BaseRNGModule() {

    vector<RNGFactory*> const &fvec = rngFactories();
    for (unsigned int i = 0; i < fvec.size(); ++i) {
	delete fvec[i];
    }
}

BaseRNGModule _base_rng_module;

