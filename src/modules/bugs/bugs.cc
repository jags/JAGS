#include <module/Module.h>
#include <distribution/DistTab.h>

#include <distributions/DBern.h>
#include <distributions/DBeta.h>
#include <distributions/DBin.h>
#include <distributions/DCat.h>
#include <distributions/DChisqr.h>
#include <distributions/DDexp.h>
#include <distributions/DDirch.h>
#include <distributions/DExp.h>
#include <distributions/DF.h>
#include <distributions/DGamma.h>
#include <distributions/DGenGamma.h>
#include <distributions/DInterval.h>
#include <distributions/DLnorm.h>
#include <distributions/DLogis.h>
#include <distributions/DMNorm.h>
#include <distributions/DMNormVC.h>
#include <distributions/DMT.h>
#include <distributions/DMulti.h>
#include <distributions/DNegBin.h>
#include <distributions/DGamPois.h>
#include <distributions/DNChisqr.h>
#include <distributions/DNorm.h>
#include <distributions/DNT.h>
#include <distributions/DPar.h>
#include <distributions/DPois.h>
#include <distributions/DRound.h>
#include <distributions/DRW1.h>
#include <distributions/DSample.h>
#include <distributions/DT.h>
#include <distributions/DUnif.h>
#include <distributions/DWeib.h>
#include <distributions/DWish.h>
#include <distributions/DHyper.h>
#include <distributions/DSum.h>
#include <distributions/SumDist.h>

#include <functions/Abs.h>
#include <functions/ArcCos.h>
#include <functions/ArcCosh.h>
#include <functions/ArcSin.h>
#include <functions/ArcSinh.h>
#include <functions/ArcTan.h>
#include <functions/ArcTanh.h>
#include <functions/BesselI.h>
#include <functions/BesselJ.h>
#include <functions/BesselK.h>
#include <functions/BesselY.h>
#include <functions/CLogLog.h>
#include <functions/Combine.h>
#include <functions/Cos.h>
#include <functions/Cosh.h>
#include <functions/DIntervalFunc.h>
#include <functions/DRoundFunc.h>
#include <functions/DSumFunc.h>
#include <functions/Exp.h>
#include <functions/Equals.h>
#include <functions/ICLogLog.h>
#include <functions/IfElse.h>
#include <functions/ILogit.h>
#include <functions/InProd.h>
#include <functions/InterpLin.h>
#include <functions/Inverse.h>
#include <functions/InverseLU.h>
#include <functions/LogDet.h>
#include <functions/LogFact.h>
#include <functions/LogGam.h>
#include <functions/Log.h>
#include <functions/Logit.h>
#include <functions/Max.h>
#include <functions/Mean.h>
#include <functions/Min.h>
#include <functions/MatMult.h>
#include <functions/Order.h>
#include <functions/Phi.h>
#include <functions/Probit.h>
#include <functions/Prod.h>
#include <functions/Rank.h>
#include <functions/Rep.h>
#include <functions/Round.h>
#include <functions/SD.h>
#include <functions/Sin.h>
#include <functions/Sinh.h>
#include <functions/Sort.h>
#include <functions/Sqrt.h>
#include <functions/Step.h>
#include <functions/Sum.h>
#include <functions/Tan.h>
#include <functions/Tanh.h>
#include <functions/Transpose.h>
#include <functions/Trunc.h>

#include <function/DFunction.h>
#include <function/PFunction.h>
#include <function/QFunction.h>

#include <samplers/CensoredFactory.h>
#include <samplers/ConjugateFactory.h>
#include <samplers/DSumFactory.h>
#include <samplers/MNormalFactory.h>
#include <samplers/DirichletFactory.h>
#include <samplers/SumFactory.h>
#include <samplers/RW1Factory.h>
#include <samplers/BinomSliceFactory.h>

using std::vector;

namespace jags {
namespace bugs {

    class BUGSModule : public Module {
    public:
	BUGSModule();
	~BUGSModule();
    };

    BUGSModule::BUGSModule() 
	: Module("bugs")
    {

	//Load distributions
	insert(new DBern);
	insert(new DBeta);
	insert(new DBin);
	insert(new DCat);
	insert(new DChisqr);
	insert(new DDexp);
	insert(new DDirch);
	insert(new DExp);
	insert(new DF);
	insert(new DGamma);
	insert(new DGenGamma);
	insert(new DLnorm);
	insert(new DLogis);
	insert(new DMNorm);
	insert(new DMNormVC);
	insert(new DMT);
	insert(new DMulti);
	insert(new DNChisqr);
	insert(new DNegBin);
	insert(new DGamPois);
	insert(new DNorm);
	insert(new DNT);
	insert(new DPar);
	insert(new DPois);
	insert(new DRW1);
	insert(new DSample);
	insert(new DT);
	insert(new DUnif);
	insert(new DWeib);
	insert(new DWish);
	insert(new DHyper);

	//Load functions
	insert(new Abs);
	insert(new ArcCos);
	insert(new ArcCosh);
	insert(new ArcSin);
	insert(new ArcSinh);
	insert(new ArcTan);
	insert(new ArcTanh);
	insert(new BesselI);
	insert(new BesselJ);
	insert(new BesselK);
	insert(new BesselY);
	insert(new Combine);
	insert(new Cos);
	insert(new Cosh);
	insert(new CLogLog);
	insert(new Exp);
	insert(new Equals);
	insert(new ICLogLog);
	insert(new IfElse);
	insert(new ILogit);
	insert(new InProd);
	insert(new InterpLin);
	insert(new Inverse);
	insert(new InverseLU);
	insert(new Log);
	insert(new LogDet);
	insert(new LogFact);
	insert(new LogGam);
	insert(new Logit);
	insert(new MatMult);
	insert(new Max);
	insert(new Mean);
	insert(new Min);
	insert(new Order);
	insert(new Phi);
	insert(new Probit);
	insert(new Prod);
	insert(new Rank);
	insert(new Rep);
	insert(new Round);
	insert(new SD);
	insert(new Sin);
	insert(new Sinh);
	insert(new Sort);
	insert(new Sqrt);
	insert(new Step);
	insert(new Sum);
	insert(new Tan);
	insert(new Tanh);
	insert(new Transpose);
	insert(new Trunc);

	//Load observable functions
	insert(new DSum, new DSumFunc);
	insert(new DInterval, new DIntervalFunc);
	insert(new DRound, new DRoundFunc);
	insert(new SumDist, new Sum);
	
	//Load sampler factories
	insert(new MNormalFactory);
	insert(new DirichletFactory);
	insert(new BinomSliceFactory);
	insert(new ConjugateFactory);
	//insert(new REFactory);
	insert(new DSumFactory);
	insert(new SumFactory);
	insert(new CensoredFactory);
	insert(new RW1Factory);
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



}}

jags::bugs::BUGSModule _bugs_module;
