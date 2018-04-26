#include <config.h>

#include "ShiftedMultinomial.h"

#include <graph/StochasticNode.h>
#include <sampler/Linear.h>
#include <sampler/SingletonGraphView.h>

#include <JRmath.h>

#include <climits>

using std::vector;
using std::string;
using std::copy;

/*
  Helper functions
*/

namespace jags {
namespace bugs {

static inline 
StochasticNode const *CHILD(SingletonGraphView const *gv, unsigned long i)
{
    return gv->stochasticChildren()[i];
}

static inline 
double SIZE(StochasticNode const *cnode, unsigned int chain)
{
    return *cnode->parents()[1]->value(chain);
}

static inline
double const *PROB(StochasticNode const *cnode, unsigned int chain)
{
    return cnode->parents()[0]->value(chain);
}

/*
  If X ~ dmulti(p, N) and X has stochastic children 
     Y[j] ~ dbin(pi[j],  n[j]) for j = 1 .. J, 
  then for each element i of X we must have either n[j] = X[i] for
  exactly one child j, or no child has a size parameter n[j] that
  depends on X.
  
  The makeIndex function checks this condition and returns a vector
  of indices such that:
     index[i] = j   if X[i] == n[j] 
     index[i] = ULONG_MAX  if no child corresponds to index i.

  If the conditions are not satisfied, an empty vector is returned,
  indicating an error.

  NB It is only safe to call this after the other checks in canSample.
*/

static vector<unsigned long>
makeIndex(SingletonGraphView const *gv, unsigned int chain)
{
    unsigned long m = gv->length(); //Length of sampled node
    vector<unsigned long> index(m, ULONG_MAX);
    unsigned long J = gv->stochasticChildren().size(); //Number of children

    // Store current value
    vector<double> xold = vector<double>(m);
    gv->getValue(xold, chain);
    
    vector<double> xnew = vector<double>(m, 0);
    gv->setValue(xnew, chain);
    for (unsigned int j = 0; j < J; ++j) {
	if (SIZE(CHILD(gv, j), chain) != 0) {
	    gv->setValue(xold, chain);
	    return vector<unsigned long>();
	}
    }

    for (unsigned long i = 0; i < m; ++i) {
	xnew[i] = 1;
	gv->setValue(xnew, chain);
	for (unsigned long j = 0; j < J; ++j) {
	    bool found_index = false;
	    double nj = SIZE(CHILD(gv,j), chain);
	    if (nj > 1) {
		// Size parameter does not correspond exactly to X[i]
		gv->setValue(xold, chain);
		return vector<unsigned long>();
	    }
	    if (nj == 1) {
		if (found_index) {
		    // More than one stochastic child depends on X[i]
		    gv->setValue(xold, chain);
		    return vector<unsigned long>();
		}
		else {
		    index[i] = j;
		    found_index = true;
		}
	    }
	}
	xnew[i] = 0;
    }

    // Restore current value
    gv->setValue(xold, chain);

    return index;
}

    ShiftedMultinomial::ShiftedMultinomial(SingletonGraphView const *gv)
	: ConjugateMethod(gv), _index(makeIndex(gv,0))
    {

    }

    bool ShiftedMultinomial::canSample(StochasticNode *snode, 
				       Graph const &graph)
    {
	/* Sampled node must have an unbounded multinomial distribution */
	if (getDist(snode) != MULTI)  
	    return false;
	if (isBounded(snode))
	    return false;

	SingletonGraphView gv(snode, graph);

	// Stochastic children must all be unbounded binomial
	vector<StochasticNode *> const &stoch_children = 
	    gv.stochasticChildren();
	for (unsigned int i = 0; i < stoch_children.size(); ++i) {
	    if (getDist(stoch_children[i]) != BIN)
		return false;
	    if (isBounded(stoch_children[i]))
		return false;
	    if (gv.isDependent(stoch_children[i]->parents()[0]))
		return false; //probability parameter depends on snode
	}
	
	// Deterministic children must be fixed scale functions
	if (!checkLinear(&gv, true))
	    return false;

	// Check indices for chain 0. See notes for MakeIndex
	vector<unsigned long> index = makeIndex(&gv, 0);
	if (index.empty()) {
	    return false; 
	}
	// Check that other chains are consistent with chain 0.
	// Possible overkill
	for (unsigned int ch = 1; ch < snode->nchain(); ++ch) {
	    if (index != makeIndex(&gv, ch))
		return false;
	}

	return true;
    }


    void ShiftedMultinomial::update(unsigned int chain, RNG *rng) const
    {
	StochasticNode *snode = _gv->node(); // Sampled node
	unsigned long m = snode->length(); // Length of sampled node
	double N = SIZE(snode, chain);
	double const *p = PROB(snode, chain);

	vector<double> prob(m);
	copy(p, p + m, prob.begin());

	vector<double> xnew = vector<double>(m, 0);
	vector<double> y = vector<double>(m, 0);

	for (unsigned long i = 0; i < m; ++i) {
	    unsigned long j = _index[i];
	    if (j != ULONG_MAX) {
		y[i] = *CHILD(_gv, j)->value(chain);
		N -= y[i];
		double pi = *PROB(CHILD(_gv, j), chain);
		prob[i] *= (1 - pi);
	    }
	}

	// Normalize probability
	double sump = 0;
	for (unsigned int i = 0; i < m; ++i) {
	    sump += prob[i];
	}
	
	// Sample multinomial as a series of binomials
	for (unsigned int i = 0; i < m - 1; i++) {
	    if (N == 0) {
		xnew[i] = 0;
	    }
	    else {
		xnew[i] = rbinom(N, prob[i]/sump, rng);
		N -= xnew[i];
		sump -= prob[i];
	    }
	}
	xnew[m - 1] = N;

	// Shift
	for (unsigned int i = 0; i < m; ++i) {
	    xnew[i] += y[i];
	}

	_gv->setValue(xnew, chain);

    }

}}
