#include <config.h>

#include <rng/RNG.h>
#include <graph/AggNode.h>
#include <graph/MixtureNode.h>
#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <sampler/Linear.h>

#include <set>
#include <stdexcept>
#include <vector>
#include <cmath>

#include "ConjugateNormal.h"
#include "ConjugateSampler.h"

#include <JRmath.h>

using std::vector;
using std::set;
using std::sqrt;
using std::invalid_argument;
using std::string;

static void calBeta(double *beta, ConjugateSampler *sampler, unsigned int chain)
{
    StochasticNode *snode = sampler->node();

    const double xold = *snode->value(chain);
    vector<StochasticNode const*> const &stoch_children = 
	sampler->stochasticChildren();

    double xnew = xold + 1;
    sampler->setValue(&xnew, 1, chain);

    double *bp = beta;    
    for (unsigned int i = 0; i < stoch_children.size(); ++i) {
	StochasticNode const *child = stoch_children[i];
	unsigned int nrow = child->length();
	double const *mu = child->parents()[0]->value(chain);
	for (unsigned int j = 0; j < nrow; ++j) {
	    bp[j] = mu[j];
	}
	bp += nrow;
    }

    sampler->setValue(&xold, 1, chain);

    bp = beta;    
    for (unsigned int i = 0; i < stoch_children.size(); ++i) {
	StochasticNode const *child = stoch_children[i];
	unsigned int nrow = child->length();
	double const *mu = child->parents()[0]->value(chain);
	for (unsigned int j = 0; j < nrow; ++j) {
	    bp[j] -= mu[j];
	}
	bp += nrow;
    }
}


ConjugateNormal::ConjugateNormal()
    : _betas(0), _length_betas(0)
{
    
}

ConjugateNormal::~ConjugateNormal()
{
    delete [] _betas;
}

void ConjugateNormal::initialize(ConjugateSampler *sampler, Graph const &graph)
{
    if (sampler->deterministicChildren().empty()) 
	return;

    vector<StochasticNode const *> const &children = 
	sampler->stochasticChildren();
    for (unsigned int i = 0; i < children.size(); ++i) {
	_length_betas += children[i]->length();
    }

    // Check for constant linear terms
    if (!checkLinear(sampler->nodes(), graph, true))
	return;

    //One-time calculation of fixed coefficients
    _betas = new double[_length_betas];
    calBeta(_betas, sampler, 0);
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
    vector<Node*> dtrm_nodes;
    Sampler::classifyChildren(vector<StochasticNode*>(1,snode), 
		              graph, stoch_nodes, dtrm_nodes);

    /* 
       Create a set of nodes containing snode and its deterministic
       descendants for the checks below.
    */
    set<Node const*> paramset;
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
	if (isBounded(stoch_nodes[i])) {
	    return false; //Truncated distribution
	}
	vector<Node const*> const &param = stoch_nodes[i]->parents();
	if (paramset.count(param[1])) {
	    return false; //Precision depends on snode
	}
    }

    // Check linearity of deterministic descendants
    if (!checkLinear(vector<StochasticNode*>(1, snode), graph, false))
	return false;

    return true; //We made it!
}

void ConjugateNormal::update(ConjugateSampler *sampler, unsigned int chain, 
			     RNG *rng) const
{
    vector<StochasticNode const*> const &stoch_children = 
	sampler->stochasticChildren();
    unsigned int nchildren = stoch_children.size();
    StochasticNode *snode = sampler->node();

    /* For convenience in the following computations, we shift the
       origin to xold, the previous value of the node */

    const double xold = *snode->value(chain);
    const double priormean = *snode->parents()[0]->value(chain) - xold; 
    const double priorprec = *snode->parents()[1]->value(chain); 

    double A = priormean * priorprec; //Weighted sum of means
    double B = priorprec; //Sum of weights

    if (sampler->deterministicChildren().empty()) {

	// This can only happen if the stochastic children are all
	// univariate normal. We know alpha = 0, beta = 1.

	for (unsigned int i = 0; i < nchildren; ++i) {
	    double Y = *stoch_children[i]->value(chain);
	    double tau = *stoch_children[i]->parents()[1]->value(chain);
	    A += (Y - xold) * tau;
	    B += tau;
	}

    }
    else {

	double *beta;
	bool temp_beta = (_betas == 0);
	if (temp_beta) {
	    beta = new double[_length_betas];
	    calBeta(beta, sampler, chain);
	}
	else {
	    beta = _betas;
	}

	double const *bp = beta;
	for (unsigned long i = 0; i < nchildren; ++i) {

	    StochasticNode const *child = stoch_children[i];

	    double const *Y = child->value(chain);
	    double const *tau = child->parents()[1]->value(chain);
	    double const *alpha = child->parents()[0]->value(chain);
	    unsigned int nrow = child->length();

	    for (unsigned int k = 0; k < nrow; ++k) {
		double tau_beta_k = 0;
		for (unsigned int k2 = 0; k2 < nrow; ++k2) {
		    tau_beta_k += tau[k * nrow + k2] * bp[k2];
		}
		A += (Y[k] - alpha[k]) * tau_beta_k;
		B += bp[k] * tau_beta_k;
	    }
	    
	    bp += nrow;
	}

	if (temp_beta) {
	    delete [] beta;
	}
    }

    // Draw the sample
    double postmean = xold + A/B;
    double postsd = sqrt(1/B);
    double xnew;

    if (isBounded(snode)) {
	Node const *lb = snode->lowerBound();
	Node const *ub = snode->upperBound();
	double plower = lb ? pnorm(*lb->value(chain), postmean, postsd, 1, 0) : 0;
	double pupper = ub ? pnorm(*ub->value(chain), postmean, postsd, 1, 0) : 1;
	double p = runif(plower, pupper, rng);
	xnew = qnorm(p, postmean, postsd, 1, 0);
    }
    else {
	xnew = rnorm(postmean, postsd, rng);  
    }
    sampler->setValue(&xnew, 1, chain);

}

string ConjugateNormal::name() const
{
    return "ConjugateNormal";
}
