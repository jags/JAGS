#include <config.h>

#include "ConjugateWishart.h"
#include "DWish.h"

#include <rng/RNG.h>
#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <graph/MixtureNode.h>
#include <sampler/Linear.h>
#include <sampler/SingletonGraphView.h>
#include <util/integer.h>

#include <set>
#include <vector>
#include <cmath>
#include <algorithm>

#include <JRmath.h>

using std::vector;
using std::set;
using std::sqrt;
using std::string;
using std::copy;

namespace jags {
namespace bugs {

static inline double getPrecision0(StochasticNode const *snode, 
				   unsigned int chain)
{
    //Returns the first element of the precision matrix for a node
    //with a multivariate normal distribution.

    return snode->parents()[1]->value(chain)[0];
}

bool ConjugateWishart::canSample(StochasticNode *snode, Graph const &graph)
{
    if (getDist(snode) != WISH)
	return false;

    if (isBounded(snode))
	return false;
  
    SingletonGraphView gv(snode, graph);
    vector<StochasticNode *> const &schild = gv.stochasticChildren();

    // Check stochastic children
    for (unsigned int i = 0; i < schild.size(); ++i) {
	if (isBounded(schild[i])) {
	    return false; //Bounded
	}
	switch(getDist(schild[i])) {
	case MNORM:
	    if (gv.isDependent(schild[i]->parents()[0])) {
		return false; //mean parameter depends on snode
	    }
	    break;
	default:
	    return false;
	}
    }

    vector<DeterministicNode*> const &dchild = gv.deterministicChildren();
    if (!dchild.empty()) {
	// Deterministic children must be scale functions
	if (!checkScale(&gv, false)) {
	    return false;
	}
	// Only mixture nodes are allowed.  If we allowed arbitrary
	// functions, the complexity would be O(nrow^4).
	for (unsigned int i = 0; i < dchild.size(); ++i) {
	    if (!isMixture(dchild[i]))
		return false;
	}
    }

    return true;
}

ConjugateWishart::ConjugateWishart(SingletonGraphView const *gv)
    : ConjugateMethod(gv)
{}

void 
ConjugateWishart::update(unsigned int chain, RNG *rng) const
{
    vector<StochasticNode *> const &stoch_children = 
	_gv->stochasticChildren();
    unsigned long nchildren = stoch_children.size();

    vector<Node const *> const &param = _gv->node()->parents();  

    double df = *param[1]->value(chain);
    double const *Rprior = param[0]->value(chain);
    unsigned long nrow = param[0]->dim()[0];

    unsigned long N = nrow * nrow;
    vector<double> R(N);
    copy(Rprior, Rprior + N, R.begin());

    //Logical mask to determine which stochastic children are active.
    vector<bool> active(nchildren, true);

    if (!_gv->deterministicChildren().empty()) {
	//Mixure model

	//Save first element of precision matrix for each child
	vector<double> precision0(nchildren); 
	for (unsigned int i = 0; i < nchildren; ++i) {
	    precision0[i] = getPrecision0(stoch_children[i], chain);
	}
	//Double the current value
	double const *x = _gv->node()->value(chain);
	vector<double> x2(N);
	for (unsigned long j = 0; j < N; ++j) {
	    x2[j] = 2 * x[j];
	}
	_gv->setValue(x2, chain);
	//See if precision matrix has changed
	for (unsigned long i = 0; i < nchildren; ++i) {
	    if (getPrecision0(stoch_children[i], chain) == precision0[i]) {
		active[i] = false; //not active
	    }
	}
    }

    for (unsigned long i = 0; i < nchildren; ++i) {
	if (active[i]) {
	    StochasticNode const *schild = stoch_children[i];
	    double const *Y = schild->value(chain);
	    double const *mu = schild->parents()[0]->value(chain);
	    
	    for (unsigned long j = 0; j < nrow; j++) {
		for (unsigned long k = 0; k < nrow; k++) {
		    R[j*nrow + k] += (Y[j] - mu[j]) * (Y[k] - mu[k]);
		}
	    }
	    df += 1;
	}
    }

    vector<double> xnew(N);

    DWish::randomSample(&xnew[0], N, &R[0], df, nrow, rng);
    _gv->setValue(xnew, chain);
}

}}
