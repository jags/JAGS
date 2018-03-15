#include <module/Module.h>

//Functions
#include <functions/Add.h>
#include <functions/And.h>
#include <functions/Divide.h>
#include <functions/Equal.h>
#include <functions/GreaterOrEqual.h>
#include <functions/GreaterThan.h>
#include <functions/LessOrEqual.h>
#include <functions/LessThan.h>
#include <functions/Multiply.h>
#include <functions/Neg.h>
#include <functions/Not.h>
#include <functions/NotEqual.h>
#include <functions/Or.h>
#include <functions/Pow.h>
#include <functions/Seq.h>
#include <functions/Subtract.h>
//Samplers
#include <samplers/SliceFactory.h>
#include <samplers/FiniteFactory.h>
//RNGs
#include <rngs/BaseRNGFactory.h>
//Monitors
#include <monitors/TraceMonitorFactory.h>
#include <monitors/MeanMonitorFactory.h>
#include <monitors/VarianceMonitorFactory.h>

using std::vector;

namespace jags {
namespace base {

    class BaseModule : public Module {

    public:
	BaseModule();
	~BaseModule();
    };

    BaseModule::BaseModule() 
	: Module("basemod")
    {

	insert(new Add);
	insert(new And);
	insert(new Divide);
	insert(new Equal);
	insert(new GreaterOrEqual);
	insert(new GreaterThan);
	insert(new LessOrEqual);
	insert(new LessThan);
	insert(new Multiply);
	insert(new Neg);
	insert(new Not);
	insert(new NotEqual);
	insert(new Or);
	insert(new Pow);
	insert(new Seq);
	insert(new Subtract);

	insert(new SliceFactory);
	insert(new FiniteFactory);
	
	insert(new BaseRNGFactory);
	
	insert(new TraceMonitorFactory);
	insert(new MeanMonitorFactory);
	insert(new VarianceMonitorFactory);
}

    BaseModule::~BaseModule() {

	vector<Function*> const &fvec = functions();
	for (unsigned int i = 0; i < fvec.size(); ++i) {
	    delete fvec[i];
	}

	vector<MonitorFactory*> const &mvec = monitorFactories();
	for (unsigned int i = 0; i < mvec.size(); ++i) {
	    delete mvec[i];
	}

	vector<RNGFactory*> const &rvec = rngFactories();
	for (unsigned int i = 0; i < rvec.size(); ++i) {
	    delete rvec[i];
	}

	vector<SamplerFactory*> const &svec = samplerFactories();
	for (unsigned int i = 0; i < svec.size(); ++i) {
	    delete svec[i];
	}
    }

}}

jags::base::BaseModule _base_module;
