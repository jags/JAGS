#include <config.h>

#include "ConjugateDirichlet.h"

#include <rng/RNG.h>
#include <graph/MixtureNode.h>
#include <graph/StochasticNode.h>

#include <set>
#include <stdexcept>
#include <vector>
#include <cmath>

#include <Rmath.h>

using std::vector;
using std::set;
using std::sqrt;
using std::invalid_argument;
using std::logic_error;

ConjugateDirichlet::ConjugateDirichlet(StochasticNode *snode, 
				       Graph const &graph,
                                       unsigned int chain)
  : ConjugateSampler(snode, graph, chain)
{
  if (!canSample(snode, graph)) {
    throw invalid_argument("Can't construct ConjugateDirichlet sampler");
  }
}

ConjugateDirichlet::~ConjugateDirichlet()
{
}

bool ConjugateDirichlet::canSample(StochasticNode *snode, Graph const &graph)
{
  if(getDist(snode) != DIRCH)
    return false;

  if (snode->isBounded())
    return false;

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

  // Check stochastic children
  for (unsigned int i = 0; i < stoch_nodes.size(); ++i) {
    vector<Node *> const &param = stoch_nodes[i]->parents();
    if (stoch_nodes[i]->isBounded()) {
      return false; //Truncated
    }
    switch(getDist(stoch_nodes[i])) {
    case CAT:
      break;
    case MULTI:
      if (param[1] == snode)
	return false;
      break;
    default:
      return false;
    }
  }

  // Check deterministic descendants
  for (unsigned int j = 0; j < dtrm_nodes.size(); ++j) {
    if (isMixture(dtrm_nodes[j])) {
      // Check that indices do not depend on snode
      if(!dtrm_nodes[j]->isLinear(paramset, false)) {
	return false;
      }
    }
    else {
      return false;
    }
  }
  
  return true;
}

static bool allzero(double const *x, long length)
{
  for (long i = 0; i < length; ++i) {
    if (x[i])
      return false;
  }
  return true;
}

void ConjugateDirichlet::update(RNG *rng)
{
    unsigned long size = node()->length();
    double *alpha = new double[size];
    double const *prior = node()->parents()[0]->value(chain());
    for (unsigned long i = 0; i < size; ++i) {
	alpha[i] = prior[i];
    }

    /* Set value to zero:
       This is an illegal value. It is used to find relevant
       children in mixture models.
       FIXME: We can ignore this in non-mixture models
    */
    double *xnew = new double[size];
    for (unsigned long i = 0; i < size; ++i) {
	xnew[i] = 0;
    }
    setValue(xnew, size);

    vector<StochasticNode const*> const &stoch_children = stochasticChildren();
    unsigned int nchildren = stoch_children.size();
    for (unsigned int i = 0; i < nchildren; ++i) {
	StochasticNode const *schild = stoch_children[i];
	long index = 0;
	double const *N = 0;
	if (allzero(schild->parents()[0]->value(chain()), 
		    schild->parents()[0]->length())) {
	    switch(_child_dist[i]) {
	    case MULTI:
	      N = schild->value(chain());
	      for (unsigned long i = 0; i < size; ++i) {
		alpha[i] += N[i];
	      }
	      break;
	    case CAT:
	      index = static_cast<long>(*schild->value(chain()) + 1.0E-6);
	      alpha[index - 1] += 1;
	      break;
	    default:
	      throw logic_error("Invalid distribution in Conjugate Dirichlet sampler");
	    }
	}
    }
    
    /* Check structural zeros */
    for (unsigned long i = 0; i < size; ++i) {
      if (prior[i] == 0 && alpha[i] != 0) {
	throw invalid_argument("Invalid likelihood for Dirichlet distribution with structural zeros");
      }
    }
  
    /* 
       Draw Dirichlet sample by drawing independent gamma random
       variates and then normalizing
    */

    double sum = 0.0;
    for (unsigned long i = 0; i < size; ++i) {
	if (alpha[i] > 0) {
	    xnew[i] = rgamma(alpha[i], 1, rng);
	    sum += xnew[i];
	}
	else {
	    xnew[i] = 0;
	}
    }
    for (unsigned long i = 0; i < size; ++i) {
	xnew[i] /= sum;
    }

    setValue(xnew, size);

    delete [] xnew;
    delete [] alpha;
}
