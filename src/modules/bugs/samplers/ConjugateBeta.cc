#include <config.h>

#include "ConjugateBeta.h"
#include "ConjugateSampler.h"

#include <rng/RNG.h>
#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <graph/MixtureNode.h>

#include <set>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <algorithm>

#include <Rmath.h>

using std::vector;
using std::set;
using std::sqrt;
using std::invalid_argument;
using std::logic_error;
using std::max;
using std::min;

ConjugateBeta::ConjugateBeta(StochasticNode *snode, Graph const &graph, unsigned int chain)
  : ConjugateSampler(snode, graph, chain) 
{
  if (!canSample(snode, graph)) {
    throw invalid_argument("Can't construct ConjugateBeta sampler");
  }
}

ConjugateBeta::~ConjugateBeta()
{

}

bool ConjugateBeta::canSample(StochasticNode *snode,
			      Graph const &graph)
{
  switch(getDist(snode)) {
  case BETA:
    break;
  case UNIF:
    // dunif(0,1) is equivalent to dbeta(1,1) 
    if(!(*snode->parents()[0]->value(0) == 0 &&
	 *snode->parents()[1]->value(0) == 1 &&
	 snode->parents()[0]->isObserved() &&
	 snode->parents()[1]->isObserved()))
      return false;
    break;
  default:
    return false;
  }

  vector<StochasticNode const*> stoch_nodes;
  vector<DeterministicNode*> dtrm_nodes;
  classifyChildren(vector<StochasticNode*>(1,snode), 
		   graph, stoch_nodes, dtrm_nodes);

  /* 
     Create a set of nodes containing snode and its deterministic
     descendants for the checks below.
  */
  set<Node *> paramset;
  paramset.insert(snode);
  for (unsigned int j = 0; j < dtrm_nodes.size(); ++j) {
    paramset.insert(dtrm_nodes[j]);
  }

  // Check deterministic descendants
  // Only Mixture nodes are allowed
  for (unsigned int j = 0; j < dtrm_nodes.size(); ++j) {
    if (isMixture(dtrm_nodes[j])) {
      // Check that indices do not depend on snode
      vector<Node*> const &parents = dtrm_nodes[j]->parents();
      unsigned int nindex = asMixture(dtrm_nodes[j])->index_size();
      for (unsigned int i = 0; i < nindex; ++i) {
        if (paramset.count(parents[i]))
          return false;
      }
    }
    else {
      return false;
    }
  }

  // Check stochastic children
  for (unsigned int i = 0; i < stoch_nodes.size(); ++i) {
    if (stoch_nodes[i]->isBounded()) {
      return false; //Bounded
    }
    switch(getDist(stoch_nodes[i])) {
    case BIN: 
      if (paramset.count(stoch_nodes[i]->parents()[1])) {
	return false; //n depends on snode
      }      
      break;
    case BERN:
      break;
    default:
      return false;
    }
  }

  return true; //We made it!
}

void ConjugateBeta::update(RNG *rng)
{
    vector<StochasticNode const*> const &stoch_children = stochasticChildren();
    StochasticNode const *snode = node();

    double a, b;
    switch (_target_dist) {
    case BETA:
	a = *snode->parents()[0]->value(chain());
	b = *snode->parents()[1]->value(chain());
	break;
    case UNIF:
	a = 1;
	b = 1;
    default:
        throw logic_error("invalid distribution in ConjugateBeta sampler");
    }
    unsigned int Nchild = stoch_children.size();

    /* For mixture models, we count only stochastic children that
       depend on snode */
    double *C = 0;
    bool is_mix = !deterministicChildren().empty();
    if (is_mix) {
	C = new double[Nchild];
	for (unsigned int i = 0; i < Nchild; ++i) {
	    C[i] = *stoch_children[i]->parents()[0]->value(chain());
	}
	// Perturb current value, keeping in the legal range [0,1]
	double x = *snode->value(chain());
	x = x > 0.5 ? x - 0.4 : x + 0.4;
	setValue(&x, 1);
	// C[i] == 1 if parameter of child i has changed (so depends on snode)
	// C[i] == 0 otherwise
	for (unsigned int i = 0; i < Nchild; ++i) {
	    C[i] = (*stoch_children[i]->parents()[0]->value(chain()) != C[i]);
	}
    }

    for (unsigned int i = 0; i < stoch_children.size(); ++i) {
	if (!(is_mix && C[i] == 0)) {
	    double y = *stoch_children[i]->value(chain());
	    double n;
	    switch(_child_dist[i]) {
	    case BIN:
		n = *stoch_children[i]->parents()[1]->value(chain());
		break;
	    case BERN:
		n = 1;
		break;
	    default:
		throw logic_error("Invalid distribution in Conjugate Beta sampler");
	    }
	    a += y;
	    b += (n - y);
	}
    }

    // Draw the sample
    double xnew = rbeta(a, b, rng);
    if (node()->isBounded()) {
	double lower = 0;
	Node const *lb = node()->lowerBound();
	if (lb) {
	    lower = max(lower, *lb->value(chain()));
	}
	double upper = 1;
	Node const *ub = node()->upperBound();
	if (ub) {
	    upper = min(upper, *ub->value(chain()));
	}
	/* Try 4 more attempts to get random sample within the bounds */
	for (int i = 0; i < 4; i++) {
	    if (xnew >= lower && xnew <= upper) {
		setValue(&xnew, 1);
		return;
	    }
	    xnew = rbeta(a, b, rng);
	}
	/* Failure! Use inversion */
	double plower = lb ? pbeta(lower, a, b, 1, 0) : 0;
	double pupper = ub ? pbeta(upper, a, b, 1, 0) : 1;
	double p = runif(plower, pupper, rng);
	xnew = qbeta(p, a, b, 1, 0);   
    }
    setValue(&xnew, 1);

    if (is_mix) {
	delete [] C;
    }

}
