#include <config.h>
#include <distribution/Distribution.h>
#include <graph/StochasticNode.h>
#include <sampler/DensitySampler.h>
#include <rng/RNG.h>
#include <util/nainf.h>
#include <graph/NodeError.h>
#include "FiniteMethod.h"

#include <cmath>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>

using std::logic_error;
using std::string;
using std::vector;
using std::exp;
using std::string;
using std::max;

namespace base {

    FiniteMethod::FiniteMethod(StochasticNode const *snode)
    {
	if (!canSample(snode)) {
	    throw logic_error("Invalid FiniteMethod");
	}

	double lower = 0, upper = 0;
	support(&lower, &upper, 1, snode, 0);

	_lower = static_cast<int>(lower);
	_upper = static_cast<int>(upper);
    }
    
    void FiniteMethod::update(RNG *rng)
    {
	int size = _upper - _lower + 1;
	vector<double> lik(size);

	//Calculate log-likelihood
	double lik_max = JAGS_NEGINF;
	for (int i = 0; i < size; i++) {
	    double ivalue = _lower + i;
	    _sampler->setValue(&ivalue, 1, _chain);
	    lik[i] = _sampler->logFullConditional(_chain);
	    lik_max = max(lik_max, lik[i]);
	}

	//Transform log-likelihood to likelihood, avoiding overflow
	double liksum = 0;
	for (int i = 0; i < size; ++i) {
	    lik[i] = exp(lik[i] - lik_max);
	    liksum += lik[i];
	}
	
	if (!jags_finite(liksum)) {
	    throw NodeError(_sampler->nodes()[0],
			    "Cannot normalize density");
	}

	/* Sample */
	double urand = rng->uniform() * liksum;
	int i;
	liksum = 0.0;
	for (i = 0; i < size - 1; i++) {
	    liksum += lik[i];
	    if (liksum > urand) {
		break;
	    }
	}
	double ivalue = _lower + i;
	_sampler->setValue(&ivalue, 1, _chain);
    }

    bool FiniteMethod::isAdaptive() const
    {
	return false;
    }

    bool FiniteMethod::adaptOff()
    {
	return true;
    }
    
    bool FiniteMethod::canSample(StochasticNode const * node)
    {
	//Node must be scalar with discrete-valued distribution of full rank
	Distribution const *dist = node->distribution();
	if (!dist->isDiscreteValued())
	    return false;

	if (node->length() != 1)
	    return false;

	if (df(node) == 0)
	    return false;

	//Support must be fixed
	if (!isSupportFixed(node)) {
	    return false;
	}
	
	//FIXME: If support is fixed it is the same for all chains
	for (unsigned int ch = 0; ch < node->nchain(); ++ch) {
	    //Distribution cannot be unbounded
	    double ulimit = JAGS_NEGINF, llimit = JAGS_POSINF;
	    support(&llimit, &ulimit, 1, node, ch);
	    if (!jags_finite(ulimit) || !jags_finite(llimit))
		return false;

	    //We don't want too many possibilities
	    double n = ulimit - llimit + 1;
	    if (n <= 1 || n > 20) //fixme: totally arbitrary
		return false;
	}
	return true;
    }

    string FiniteMethod::name() const
    {
	return "FiniteMethod";
    }
}
