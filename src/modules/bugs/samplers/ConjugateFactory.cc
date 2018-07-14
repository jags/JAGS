#include <config.h>

#include "ConjugateFactory.h"
#include "ConjugateNormal.h"
#include "ConjugateGamma.h"
#include "ConjugateBeta.h"
#include "ConjugateDirichlet.h"
#include "ConjugateMNormal.h"
#include "ConjugateWishart.h"
//#include "TruncatedGamma.h"
#include "ShiftedCount.h"
#include "ShiftedMultinomial.h"

#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>
#include <sampler/SingletonGraphView.h>
#include <module/ModuleError.h>
#include <sampler/ImmutableSampler.h>

#include <string>

using std::string;

namespace jags {
namespace bugs {

bool ConjugateFactory::canSample(StochasticNode * snode,
				 Graph const &graph) const
{
    bool ans = false;
    switch(getDist(snode)) {
    case NORM:
	ans = ConjugateNormal::canSample(snode, graph);
	break;
    case GAMMA: case CHISQ:
	ans = ConjugateGamma::canSample(snode, graph);
	break;
    case EXP:
	ans = ConjugateGamma::canSample(snode, graph) ||
	    ConjugateNormal::canSample(snode, graph);
	break;
    case BETA:
	ans = ConjugateBeta::canSample(snode, graph);
	break;
    case DIRCH:
	ans = ConjugateDirichlet::canSample(snode, graph);
	break;
    case MNORM:
	ans = ConjugateMNormal::canSample(snode, graph);
	break;
    case WISH:
	ans = ConjugateWishart::canSample(snode, graph);
	break;
    case UNIF:
        /*
	ans = TruncatedGamma::canSample(snode, graph) ||
	      ConjugateBeta::canSample(snode, graph);
        */
	ans = ConjugateBeta::canSample(snode, graph);
	break;
    case POIS: case BIN: case NEGBIN:
	ans = ShiftedCount::canSample(snode, graph);
	break;
    case MULTI:
	ans = ShiftedMultinomial::canSample(snode, graph);
	break;
    case BERN: case CAT: case DEXP: case LNORM: case LOGIS:  case PAR: 
    case T: case WEIB: case OTHERDIST:
	break;
    }
    
    return ans;
}

Sampler *ConjugateFactory::makeSampler(StochasticNode *snode, 
				       Graph const &graph) const
{
    SingletonGraphView *gv = new SingletonGraphView(snode, graph);
    ConjugateMethod* method = nullptr;
    string name;

    switch (getDist(snode)) {
    case NORM:
	method = new ConjugateNormal(gv);
	name = "bugs::ConjugateNormal";
	break;
    case GAMMA: case CHISQ:
	method = new ConjugateGamma(gv);
	name = "bugs::ConjugateGamma";
	break;
    case EXP:
	if (ConjugateGamma::canSample(snode, graph)) {
	    method = new ConjugateGamma(gv);
	    name = "bugs::ConjugateGamma";
	}
	else if (ConjugateNormal::canSample(snode, graph)) {
	    method = new ConjugateNormal(gv);
	    name = "bugs::ConjugateNormal";
	}
	else {
	    throwLogicError("Cannot find conjugate sampler for exponential");
	}
	break;
    case BETA:
	method = new ConjugateBeta(gv);
	name = "bugs::ConjugateBeta";
	break;
    case DIRCH:
	method = new ConjugateDirichlet(gv);
	name = "bugs::ConjugateDirichlet";
	break;
    case MNORM:
	method = new ConjugateMNormal(gv);
	name = "bugs::ConjugateMNormal";
	break;
    case WISH:
	method = new ConjugateWishart(gv);
	name = "bugs::ConjugateWishart";
	break;
    case UNIF:
	/*
	  if (TruncatedGamma::canSample(snode, graph)) {
	  method = new TruncatedGamma(gv);
	  name = "TruncatedGamma";
	  }
	  else 
	*/
	if (ConjugateBeta::canSample(snode, graph)) {
	    method = new ConjugateBeta(gv);
	    name = "bugs::ConjugateBeta";
	}
	else {
	    throwLogicError("Cannot find conjugate sampler for uniform");
	}
	break;
    case POIS: case BIN: case NEGBIN:
	method = new ShiftedCount(gv);
	name = "bugs::ShiftedCount";
	break;
    case MULTI:
	method = new ShiftedMultinomial(gv);
	name = "bugs::ShiftedMultiNomial";
	break;
    case BERN: case CAT: case DEXP: case LNORM: case LOGIS:  case PAR: 
    case T: case WEIB: case OTHERDIST:
	throwLogicError("Unable to create conjugate sampler");
    }
    
    return new ImmutableSampler(gv, method, name);
}

string ConjugateFactory::name() const
{
    return "bugs::Conjugate";
}

}}
