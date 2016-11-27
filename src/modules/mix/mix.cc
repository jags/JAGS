#include <module/Module.h>

#include <distributions/DNormMix.h>
#include <distributions/DBetaBin.h>
#include <distributions/DPick.h>
#include <samplers/MixSamplerFactory.h>
#include <samplers/DirichletCatFactory.h>
#include <samplers/LDAFactory.h>

using std::vector;

namespace jags {
namespace mix {

    class MIXModule : public Module {
    public:
	MIXModule();
	~MIXModule();
    };

    MIXModule::MIXModule() 
	: Module("mix") {

	insert(new DBetaBin);
	insert(new DNormMix);
        insert(new DPick);
	
	insert(new MixSamplerFactory);
	insert(new DirichletCatFactory);
	insert(new LDAFactory);
    }

    MIXModule::~MIXModule() {

	vector<Distribution*> const &dvec = distributions();
	for (unsigned int i = 0; i < dvec.size(); ++i) {
	    delete dvec[i];
	}

	vector<SamplerFactory*> const &svec = samplerFactories();
	for (unsigned int i = 0; i < svec.size(); ++i) {
	    delete svec[i];
	}
    }

}}

jags::mix::MIXModule _mix_module;
