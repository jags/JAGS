#include <config.h>

#include "ConjugateBeta.h"

#include <rng/RNG.h>
#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <graph/MixtureNode.h>
#include <sampler/SingletonGraphView.h>
#include <sampler/Linear.h>
#include <module/ModuleError.h>
#include <DBeta.h>

#include <set>
#include <vector>
#include <cmath>
#include <algorithm>

#include <JRmath.h>

#include <iostream>

using std::vector;
using std::set;
using std::sqrt;
using std::max;
using std::min;
using std::string;

namespace jags {
namespace bugs {

bool ConjugateBeta::canSample(StochasticNode *snode, Graph const &graph)
{
    ConjugateDist dist = getDist(snode);
    if (dist !=  UNIF && dist != BETA) {
	return false;
    }

    SingletonGraphView gv(snode, graph);
    vector<DeterministicNode*> const &dchild = gv.deterministicChildren();
    vector<StochasticNode *> const &schild = gv.stochasticChildren();

    // Check deterministic descendants
    // Only Mixture nodes are allowed
    for (unsigned int j = 0; j < dchild.size(); ++j) {
	if (!isMixture(dchild[j])) {
	    return false;
	}
    }
    if (!checkScale(&gv, false)) {
	return false;
    }

    // Check stochastic children
    for (unsigned long i = 0; i < schild.size(); ++i) {
	if (isBounded(schild[i])) {
	    return false; //Bounded child
	}
	switch(getDist(schild[i])) {
	case BIN: case NEGBIN:
	    if (gv.isDependent(schild[i]->parents()[1])) {
		return false; //size parameter depends on snode
	    }      
	    break;
	case BERN:
	    break;
	case BETA: case CAT: case CHISQ: case DEXP: case DIRCH: case EXP:
	case GAMMA: case LNORM: case LOGIS: case MNORM: case MULTI: case NORM:
	case PAR: case POIS: case T: case UNIF: case WEIB: case WISH:
	case OTHERDIST:
	    return false;
	}
    }

    return true; //We made it!
}


ConjugateBeta::ConjugateBeta(SingletonGraphView const *gv)
    : ConjugateMethod(gv)
{
}

void ConjugateBeta::update(unsigned int chain, RNG *rng) const
{
    vector<StochasticNode *> const &stoch_children = 
	_gv->stochasticChildren();
    StochasticNode const *snode = _gv->node();

    double a=0, b=0; //-Wall
    switch (_target_dist) {
    case BETA:
	a = *snode->parents()[0]->value(chain);
	b = *snode->parents()[1]->value(chain);
	break;
    case UNIF:
	a = 1;
	b = 1;
	break;
    case BERN: case BIN: case CAT: case CHISQ: case DEXP: case DIRCH: case EXP:
    case GAMMA: case LNORM: case LOGIS: case MNORM: case MULTI: case NEGBIN:
    case NORM: case PAR: case POIS: case T: case WEIB: case WISH:
    case OTHERDIST:
	throwLogicError("Invalid distribution in ConjugateBeta sampler");
    }
    unsigned long Nchild = stoch_children.size();

    /* Get bounds */
    double lower = 0, upper = 1;
    snode->support(&lower, &upper, 1, chain);
    /* Clamp bounds to (0,1) */
    lower = max(lower, 0.0);
    upper = min(upper, 1.0);
    if (lower >= upper) {
	throwNodeError(snode, "Invalid bounds in ConjugateBeta");
    }
    
    /* For mixture models, we count only stochastic children that
       depend on snode */
    bool is_mix = !_gv->deterministicChildren().empty();
    vector<double> C(Nchild, 1);
    if (is_mix) {
	for (unsigned long i = 0; i < Nchild; ++i) {
	    C[i] = *stoch_children[i]->parents()[0]->value(chain);
	}
	// Perturb current value, keeping in the legal range
	double midpoint = (lower + upper)/2;
	double x = *snode->value(chain);
	x = x > midpoint ? x - (x - lower)/2 : x + (upper - x)/2;
	_gv->setValue(&x, 1, chain);
	// C[i] == 1 if parameter of child i has changed (so depends on snode)
	// C[i] == 0 otherwise
	for (unsigned long i = 0; i < Nchild; ++i) {
	    C[i] = (*stoch_children[i]->parents()[0]->value(chain) != C[i]);
	}
    }

    for (unsigned long i = 0; i < stoch_children.size(); ++i) {
	if (!(is_mix && C[i] == 0)) {
	    double y = *stoch_children[i]->value(chain);
	    double n;
	    switch(_child_dist[i]) {
	    case BIN:
		n = *stoch_children[i]->parents()[1]->value(chain);
		a += y;
		b += n - y;
		break;
	    case NEGBIN:
		n = *stoch_children[i]->parents()[1]->value(chain);
		a += n;
		b += y;
		break;
	    case BERN:
		a += y;
		b += 1 - y;
		break;
	    case BETA: case CAT: case CHISQ: case DEXP: case DIRCH: case EXP:
	    case GAMMA: case LNORM: case LOGIS: case MNORM: case MULTI:
	    case NORM: case PAR: case POIS: case T: case UNIF: case WEIB:
	    case WISH: case OTHERDIST:
		throwLogicError("Invalid distribution in Conjugate Beta sampler");
	    }
	}
    }

    // Draw the sample
    DBeta dbeta;
    vector<double const *> par = { &a, &b };
    double xnew = dbeta.randomSample(par, &lower, &upper, rng);
    _gv->setValue(&xnew, 1, chain);
}

}}
