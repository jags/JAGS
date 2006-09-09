#include <config.h>

#include <rng/RNG.h>
#include <graph/AggNode.h>
#include <graph/MixtureNode.h>
#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>

#include <set>
#include <stdexcept>
#include <vector>
#include <cmath>

#include "ConjugateNormal.h"
#include "ConjugateSampler.h"

#include <Rmath.h>

using std::vector;
using std::set;
using std::sqrt;
using std::invalid_argument;

ConjugateNormal::ConjugateNormal(StochasticNode *snode, Graph const &graph, 
				 unsigned int chain)
    : ConjugateSampler(snode, graph, chain), _betas(0), _length_betas(0)
{
    if (!canSample(snode, graph)) {
	throw invalid_argument("Can't construct ConjugateNormal sampler");
    }

    if (!deterministicChildren().empty()) {
	vector<StochasticNode const *> const &children = stochasticChildren();
	for (unsigned int i = 0; i < children.size(); ++i) {
	    _length_betas += children[i]->length();
	}
    
	// Check for constant linear terms
	set<Node*> paramset;
	vector<DeterministicNode*> const &dtrm = deterministicChildren();
	paramset.insert(snode);
	for (unsigned int j = 0; j < dtrm.size(); ++j) {
	    paramset.insert(dtrm[j]);
	}

	bool fixed_beta = true;
	for (unsigned int j = 0; j < dtrm.size(); ++j) {
	    if (!dtrm[j]->isLinear(paramset, true)) {
		fixed_beta = false;
		break;
	    }
	}
	
	if (fixed_beta) {
	    _betas = new double[_length_betas];
	    calBeta();
	}
    }
}

ConjugateNormal::~ConjugateNormal()
{
    delete [] _betas;
}

void ConjugateNormal::calBeta()
{
    const double xold = *node()->value(chain());
    vector<StochasticNode const*> const &stoch_children = stochasticChildren();

    double xnew = xold + 1;
    setValue(&xnew, 1);

    double *beta = _betas;    
    for (unsigned int i = 0; i < stoch_children.size(); ++i) {
	StochasticNode const *snode = stoch_children[i];
	unsigned int nrow = snode->length();
	double const *mu = snode->parents()[0]->value(chain());
	for (unsigned int j = 0; j < nrow; ++j) {
	    beta[j] = mu[j];
	}
	beta += nrow;
    }

    setValue(&xold, 1);

    beta = _betas;    
    for (unsigned int i = 0; i < stoch_children.size(); ++i) {
	StochasticNode const *snode = stoch_children[i];
	unsigned int nrow = snode->length();
	double const *mu = snode->parents()[0]->value(chain());
	for (unsigned int j = 0; j < nrow; ++j) {
	    beta[j] -= mu[j];
	}
	beta += nrow;
    }
}


bool ConjugateNormal::canSample(StochasticNode *snode, Graph const &graph)
{
  /*
    1) Stochastic children of sampled node must be normal, must not
    be truncated, and must depend on snode only via the mean parameter
    2) The mean parameter must be a linear function of snode. 
  */

  if (getDist(snode) != NORM)
    return false;

  vector<StochasticNode const*> stoch_nodes;
  vector<DeterministicNode*> dtrm_nodes;
  classifyChildren(vector<StochasticNode*>(1,snode), 
		   graph, stoch_nodes, dtrm_nodes);

  /* 
     Create a set of nodes containing snode and its deterministic
     descendants for the checks below.
  */
  set<Node*> paramset;
  paramset.insert(snode);
  for (unsigned int j = 0; j < dtrm_nodes.size(); ++j) {
    paramset.insert(dtrm_nodes[j]);
  }

  // Check stochastic children
  for (unsigned int i = 0; i < stoch_nodes.size(); ++i) {
      switch (getDist(stoch_nodes[i])) {
      case NORM: case MNORM:
	  break;
      default:
	  return false; //Not normal
      }
      if (stoch_nodes[i]->isBounded()) {
	  return false; //Truncated distribution
      }
      vector<Node*> const &param = stoch_nodes[i]->parents();
      if (paramset.count(param[1])) {
	  return false; //Precision depends on snode
      }
  }
  
  // Check linearity of deterministic descendants
  for (unsigned int j = 0; j < dtrm_nodes.size(); ++j) {
    if (!dtrm_nodes[j]->isLinear(paramset, false))
      return false;
  }

  return true; //We made it!
}


void ConjugateNormal::update(RNG *rng)
{
    vector<StochasticNode const*> const &stoch_children = stochasticChildren();
    unsigned int nchildren = stoch_children.size();

    /* For convenience in the following computations, we shift the
       origin to xold, the previous value of the node */

    const double xold = *node()->value(chain());
    const double priormean = *node()->parents()[0]->value(chain()) - xold; 
    const double priorprec = *node()->parents()[1]->value(chain()); 

    double A = priormean * priorprec; //Weighted sum of means
    double B = priorprec; //Sum of weights

    if (deterministicChildren().empty()) {

	// This can only happen if the stochastic children are all
	// univariate normal. We know alpha = 0, beta = 1.

	for (unsigned int i = 0; i < nchildren; ++i) {
	    double Y = *stoch_children[i]->value(chain());
	    double tau = *stoch_children[i]->parents()[1]->value(chain());
	    A += (Y - xold) * tau;
	    B += tau;
	}

    }
    else {

	bool temp_beta = (_betas == 0);
	if (temp_beta) {
	    _betas = new double[_length_betas];
	    calBeta();
	}

	double const *beta = _betas;
	for (unsigned long i = 0; i < nchildren; ++i) {

	    StochasticNode const *snode = stoch_children[i];

	    double const *Y = snode->value(chain());
	    double const *tau = snode->parents()[1]->value(chain());
	    double const *alpha = snode->parents()[0]->value(chain());
	    unsigned int nrow = snode->length();

	    for (unsigned int k = 0; k < nrow; ++k) {
		double tau_beta_k = 0;
		for (unsigned int k2 = 0; k2 < nrow; ++k2) {
		    tau_beta_k += tau[k * nrow + k2] * beta[k2];
		}
		A += (Y[k] - alpha[k]) * tau_beta_k;
		B += beta[k] * tau_beta_k;
	    }
	    
	    beta += nrow;
	}

	if (temp_beta) {
	    delete [] _betas;
	    _betas = 0;
	}
    
    }

    // Draw the sample
    double postmean = xold + A/B;
    double postsd = sqrt(1/B);
    double xnew;

    if (node()->isBounded()) {
	Node const *lb = node()->lowerBound();
	Node const *ub = node()->upperBound();
	double plower = lb ? pnorm(*lb->value(chain()), postmean, postsd, 1, 0) : 0;
	double pupper = ub ? pnorm(*ub->value(chain()), postmean, postsd, 1, 0) : 1;
	double p = runif(plower, pupper, rng);
	xnew = qnorm(p, postmean, postsd, 1, 0);
    }
    else {
	xnew = rnorm(postmean, postsd, rng);  
    }
    setValue(&xnew, 1);

}
