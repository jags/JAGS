#include "OrderedProbit.h"
#include "Classify.h"

#include "KS.h"

#include <graph/StochasticNode.h>
#include <rng/TruncatedNormal.h>
#include <rng/RNG.h>

#include <cmath>

using std::exp;
using std::log;
using std::sqrt;

namespace jags {
namespace glm {

    OrderedProbit::OrderedProbit(StochasticNode const *snode,
				 unsigned int chain)
	: Outcome(snode, chain), _y(snode->value(chain)[0]),
	  _cuts(snode->parents()[1]->value(chain)),
	  _ncut(snode->parents()[1]->length()), _z(0)
    {
	//fixme: sanity checks on snode
    }

    double OrderedProbit::value() const 
    {
	return _z;
    }

    double OrderedProbit::precision() const 
    {
	return 1;
    }
    
    void OrderedProbit::update(RNG *rng)
    {
	/* Albert Chib update */
	int x = static_cast<int>(_y) - 1;
	if (x == 0) {
	    _z = rnormal(_cuts[0], rng, _lp);
	}
	else if (x == _ncut) {
	    _z = lnormal(_cuts[_ncut - 1], rng, _lp);
	}
	else {
	    _z = inormal(_cuts[x-1], _cuts[x], rng, _lp);
	}
    }

    void OrderedProbit::update(double mean, double var, RNG *rng)
    {
	/* Holmes-Held update */
	int x = static_cast<int>(_y) - 1;
	if (x == 0) {
	    _z = rnormal(_cuts[0], rng, mean, sqrt(var + 1));
	}
	else if (x == _ncut) {
	    _z = lnormal(_cuts[_ncut - 1], rng, mean, sqrt(var + 1));
	}
	else {
	    _z = inormal(_cuts[x-1], _cuts[x], rng, mean, sqrt(var + 1));
	}
    }

    bool OrderedProbit::canRepresent(StochasticNode const *snode)
    {
	return getFamily(snode) == GLM_ORDPROBIT &&
	    getLink(snode) == LNK_LINEAR;
    }

}}
