#include <config.h>

#include "ShiftedCount.h"

#include <graph/StochasticNode.h>
#include <sampler/SingletonGraphView.h>
#include <module/ModuleError.h>

#include <JRmath.h>

using std::vector;
using std::string;

namespace jags {
namespace bugs {

    ShiftedCount::ShiftedCount(SingletonGraphView const *gv)
	: ConjugateMethod(gv)
    {
    }

    bool ShiftedCount::canSample(StochasticNode *snode, Graph const &graph)
    {
	switch(getDist(snode)) {
	case POIS: case BIN: case NEGBIN:
	    break;
	case BERN: case BETA: case CAT: case CHISQ: case DEXP: case DIRCH:
	case EXP: case GAMMA: case LNORM: case LOGIS: case MNORM: case MULTI:
	case NORM: case PAR: case T: case UNIF: case WEIB: case WISH:
	case OTHERDIST:
	    return false;
	}
	if (isBounded(snode))
	    return false;

	SingletonGraphView gv(snode, graph);

	// There should be no deterministic children
	if (!gv.deterministicChildren().empty())
	    return false;

	// Check stochastic children: Must have single child node ...
	vector<StochasticNode *> const &child_nodes = 
	    gv.stochasticChildren();
	if (child_nodes.size() != 1) 
	    return false;
	
	// ... with a Binomial(p, N) distribution such that N is the
	// sampled node and p does not depend on sampled node
	StochasticNode const *cnode = child_nodes.front();
	if (getDist(cnode) != BIN) 
	    return false;
	if (isBounded(cnode))
	    return false;
	if (cnode->parents()[1] != snode)
	    return false;
	if (cnode->parents()[0] == snode)
	    return false;

	return true;
    }


    void ShiftedCount::update(unsigned int chain, RNG *rng) const
    {
	StochasticNode const* snode = _gv->node();
	StochasticNode * cnode = _gv->stochasticChildren().front();

	double y = *cnode->value(chain); //Child value
	// Prior rate (Poisson) or probability (Binomial)
	double lambda = *snode->parents()[0]->value(chain); 
	// Binomial probability of child
	double pi = *cnode->parents()[0]->value(chain); 

	double xnew = y;
	double N = 0;
	switch(_target_dist) {
	case POIS:
	    lambda = (1 - pi) * lambda;
	    xnew += rpois(lambda, rng);
	    break;
	case BIN:
	    N = *snode->parents()[1]->value(chain);
	    lambda = (1 - pi) * lambda  / ((1 - pi) * lambda + (1 - lambda));
	    xnew += rbinom(N - y, lambda, rng);
	    break;
	case NEGBIN:
	    N = *snode->parents()[1]->value(chain);
	    lambda = (1 - pi) * lambda + pi;
	    xnew += rnbinom(N - y, lambda, rng);
	    break;
	case BERN: case BETA: case CAT: case CHISQ: case DEXP: case DIRCH:
	case EXP: case GAMMA: case LNORM: case LOGIS: case MNORM: case MULTI:
	case NORM: case PAR: case T: case UNIF: case WEIB: case WISH:
	case OTHERDIST:
	    throwLogicError("Invalid distribution in ShiftedCount");
	}
	_gv->setValue(&xnew, 1, chain);
    }

}}
