//#define bugs_LTX_load load

#include <Module.h>
#include <distribution/DistTab.h>

#include <distributions/DBern.h>
#include <distributions/DBeta.h>
#include <distributions/DBin.h>
#include <distributions/DCat.h>
#include <distributions/DChisqr.h>
#include <distributions/DDexp.h>
#include <distributions/DDirch.h>
#include <distributions/DExp.h>
#include <distributions/DGamma.h>
#include <distributions/DInterval.h>
#include <distributions/DLnorm.h>
#include <distributions/DLogis.h>
#include <distributions/DMNorm.h>
#include <distributions/DMState.h>
#include <distributions/DMulti.h>
#include <distributions/DNegBin.h>
#include <distributions/DNorm.h>
#include <distributions/DPar.h>
#include <distributions/DPois.h>
#include <distributions/DT.h>
#include <distributions/DUnif.h>
#include <distributions/DWeib.h>
#include <distributions/DWish.h>
#include <distributions/DHyper.h>
#include <distributions/DSum.h>
#include <distributions/DNormMix.h>

#include <functions/CLogLog.h>
#include <functions/Exp.h>
#include <functions/ICLogLog.h>
#include <functions/ILogit.h>
#include <functions/InProd.h>
#include <functions/Inverse.h>
#include <functions/LogDet.h>
#include <functions/LogFact.h>
#include <functions/LogGam.h>
#include <functions/Log.h>
#include <functions/Logit.h>
#include <functions/Max.h>
#include <functions/Mean.h>
#include <functions/Mexp.h>
#include <functions/Min.h>
#include <functions/Phi.h>
#include <functions/Probit.h>
#include <functions/Rank.h>
#include <functions/SD.h>
#include <functions/Sort.h>
#include <functions/Sqrt.h>
#include <functions/Step.h>
#include <functions/Sum.h>

#include <samplers/ConjugateFactory.h>
#include <samplers/DSumFactory.h>

using std::vector;

class BUGSModule : public Module {
public:
    BUGSModule();
    ~BUGSModule();
};

BUGSModule::BUGSModule() {

    //Load distributions
    insert(new DBern);
    insert(new DBeta);
    insert(new DBin);
    insert(new DCat);
    insert(new DChisqr);
    insert(new DDexp);
    insert(new DDirch);
    insert(new DExp);
    insert(new DGamma);
    insert(new DInterval);
    insert(new DLnorm);
    insert(new DLogis);
    insert(new DMNorm);
    insert(new DMState);
    insert(new DMulti);
    insert(new DNegBin);
    insert(new DNorm);
    insert(new DPar);
    insert(new DPois);
    insert(new DT);
    insert(new DUnif);
    insert(new DWeib);
    insert(new DWish);
    insert(new DHyper);
    insert(new DSum);
    insert(new DNormMix);

    insert(new CLogLog);
    insert(new Exp);
    insert(new ICLogLog);
    insert(new ILogit);
    insert(new InProd);
    insert(new Inverse);
    insert(new Log);
    insert(new LogDet);
    insert(new LogFact);
    insert(new LogGam);
    insert(new Logit);
    insert(new Max);
    insert(new Mean);
    insert(new Mexp);
    insert(new Min);
    insert(new Phi);
    insert(new Probit);
    insert(new Rank);
    insert(new SD);
    insert(new Sort);
    insert(new Sqrt);
    insert(new Step);
    insert(new Sum);

    insert(new ConjugateFactory);
    insert(new DSumFactory);
}

BUGSModule::~BUGSModule() {

    vector<Distribution*> const &dvec = distributions();
    for (unsigned int i = 0; i < dvec.size(); ++i) {
	delete dvec[i];
    }

    vector<Function*> const &fvec = functions();
    for (unsigned int i = 0; i < fvec.size(); ++i) {
	delete fvec[i];
    }

    vector<SamplerFactory*> const &svec = samplerFactories();
    for (unsigned int i = 0; i < svec.size(); ++i) {
	delete svec[i];
    }
}

/*
extern "C" {
    void load() {
	static BUGSModule bugs_module;
	return;
    }
}
*/

BUGSModule _bugs_module;
