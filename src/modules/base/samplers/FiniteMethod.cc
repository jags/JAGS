#include <config.h>
#include <distribution/Distribution.h>
#include <graph/StochasticNode.h>
#include <rng/RNG.h>
#include <util/nainf.h>
#include <graph/NodeError.h>
#include <sampler/SingletonGraphView.h>
#include <module/ModuleError.h>

#include "FiniteMethod.h"

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

using std::string;
using std::vector;
using std::exp;
using std::string;
using std::binary_search;
using std::upper_bound;

static double mkLower(jags::SingletonGraphView const *gv) {
    double lower = 0, upper = 0;
    gv->node()->support(&lower, &upper, 1, 0);
    return lower;
}

static double mkUpper(jags::SingletonGraphView const *gv) {
    double lower = 0, upper = 0;
    gv->node()->support(&lower, &upper, 1, 0);
    return upper;
}

static unsigned long mkSize(jags::SingletonGraphView const *gv) {
    double lower = 0, upper = 0;
    gv->node()->support(&lower, &upper, 1, 0);
    return static_cast<unsigned long>(upper - lower + 1);
}

namespace jags {
namespace base {

    FiniteMethod::FiniteMethod(SingletonGraphView const *gv)
	: _gv(gv), _lower(mkLower(gv)), _upper(mkUpper(gv)), _size(mkSize(gv))
    {
	if (!canSample(gv->node())) {
	    throwLogicError("Invalid FiniteMethod");
	}
    }
    
    void FiniteMethod::update(unsigned int chain, RNG *rng) const
    {
	vector<double> lik(_size);

	//Calculate log-likelihood
	double lik_max = JAGS_NEGINF;
	for (unsigned long i = 0; i < _size; i++) {
	    double ivalue = _lower + i;
	    _gv->setValue(&ivalue, 1, chain);
	    lik[i] = _gv->logFullConditional(chain);
	    if (lik[i] > lik_max) lik_max = lik[i];
	}
	
	//Transform log-likelihood to likelihood, avoiding overflow
	//and calculate partial sums
	double liksum = 0;
	for (unsigned long i = 0; i < _size; ++i) {
	    liksum += exp(lik[i] - lik_max);
	    lik[i] = liksum;
	}
	
	if (!jags_finite(liksum)) {
	    throwNodeError(_gv->node(), "Cannot normalize density");
	}

	// Sample
	double urand = rng->uniform() * liksum;
	long i = upper_bound(lik.begin(), lik.end(), urand) - lik.begin();
	double ivalue = _lower + i;
	_gv->setValue(&ivalue, 1, chain);
    }

    bool FiniteMethod::canSample(StochasticNode const * node)
    {
	//Node must be scalar with discrete-valued distribution of full rank
	if (!node->isDiscreteValued())
	    return false;

	if (node->length() != 1)
	    return false;

	if (!node->fullRank())
	    return false;

	//Support must be fixed
	if (!isSupportFixed(node)) {
	    return false;
	}

	//Distribution cannot be unbounded
	double ulimit = JAGS_NEGINF, llimit = JAGS_POSINF;
	node->support(&llimit, &ulimit, 1, 0);
	if (!jags_finite(ulimit) || !jags_finite(llimit))
	    return false;
	
	/* If there are too many possibilities then FiniteMethod,
	   becomes slow as it has to evaluate the likelihood at each
	   point (inversion). So we drop FiniteMethod in favour of the
	   discrete slice sampler if the number of possibilities is
	   too large. It's not clear what "too large" is. Here we use
	   a cutoff size of 100, unless the node has a categorical
	   prior, noting that dcat can take an arbitrary, unordered
	   vector of probabilities as argument.
	   
	   FIXME: a more careful criterion would be to test whether
	   the log-density is unimodal.
	*/
	if (node->distribution()->name() != "dcat") {
	    if (ulimit - llimit >= 100)
		return false;
	}

	return true;
    }

}}
