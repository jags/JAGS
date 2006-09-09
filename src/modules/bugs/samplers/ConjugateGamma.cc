#include <config.h>

#include "ConjugateGamma.h"

#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <graph/MixtureNode.h>
#include <graph/NodeError.h>
#include <sarray/SArray.h>

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
using std::sort;

ConjugateGamma::ConjugateGamma(StochasticNode *snode, Graph const &graph, unsigned int chain)
    : ConjugateSampler(snode, graph, chain), _coef(0)
{
    if (!canSample(snode, graph)) {
	throw invalid_argument("Can't construct ConjugateGamma sampler");
    }

    if (!deterministicChildren().empty()) {
    
	// Check for constant scale transformation
	set<Node*> paramset;
	vector<DeterministicNode*> const &dtrm = deterministicChildren();
	paramset.insert(snode);
	for (unsigned int j = 0; j < dtrm.size(); ++j) {
	    paramset.insert(dtrm[j]);
	}

	bool fixed_coef = true;
	for (unsigned int j = 0; j < dtrm.size(); ++j) {
	    if (!dtrm[j]->isScale(paramset, true)) {
		fixed_coef = false;
		break;
	    }
	}
	
	if (fixed_coef) {
	    _coef = new double[stochasticChildren().size()];
	    calCoef();
	}
    }
}

ConjugateGamma::~ConjugateGamma()
{
    delete [] _coef;
}

static double 
getScale(StochasticNode const *snode, ConjugateDist d, unsigned int chain)
{
    //Get scale parameter of snode
    switch(d) {
    case GAMMA: case NORM: case DEXP: case WEIB:
	return *snode->parents()[1]->value(chain);
	break;
    case EXP: case POIS:
	return *snode->parents()[0]->value(chain);
	break;
    default:
	throw NodeError(snode, 
			"Can't get scale parameter: invalid distribution");
    } 
}

void ConjugateGamma::calCoef()
{
    const double xold = *node()->value(chain());
    vector<StochasticNode const*> const &stoch_children = stochasticChildren();
    unsigned long nchildren = stoch_children.size();

    for (unsigned int i = 0; i < nchildren; ++i) {
	_coef[i] = -getScale(stoch_children[i], _child_dist[i], chain());
    }
    double val = xold + 1;
    setValue(&val, 1);
    for (unsigned int i = 0; i < nchildren; ++i) {
	_coef[i] += getScale(stoch_children[i], _child_dist[i], chain());
    }
    setValue(&xold, 1);
}


bool ConjugateGamma::canSample(StochasticNode *snode, Graph const &graph)
{
    switch (getDist(snode)) {
    case GAMMA: case EXP: case CHISQ:
	/* 
	   The exponential and chisquare distributions are both special
	   cases of the gamma distribution and are handled by the conjugate
	   gamma sampler.
	*/
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

    // Check stochastic children
    for (unsigned int i = 0; i < stoch_nodes.size(); ++i) {
	vector<Node*> const &param = stoch_nodes[i]->parents();
	if (stoch_nodes[i]->isBounded()) {
	    return false; //Bounded
	}
	switch(getDist(stoch_nodes[i])) {
	case EXP: case POIS:
	    break;
	case GAMMA: case NORM: case DEXP: case WEIB:
	    if (paramset.count(param[0])) {
		return false; //non-scale parameter depends on snode
	    }
	    break;
	default:
	    return false;
	}
    }

  // Check deterministic descendants are scale transformations 
  for (unsigned int j = 0; j < dtrm_nodes.size(); ++j) {
    if (!dtrm_nodes[j]->isScale(paramset, false))
      return false;
  }

    /*
    // Check deterministic descendants
    for (unsigned int j = 0; j < dtrm_nodes.size(); ++j) {
	LogicalNode const *lnode = asLogical(dtrm_nodes[j]);
	if (lnode) {
	    unsigned int nfactor = 0;
	    vector<Node*> const &param = lnode->parents();
	    switch(getOp(lnode)) {
	    case DIVIDE:
		if (paramset.count(param[1]))
		    return false; //reciprocal term
		break;
	    case MULTIPLY:
		for (unsigned int k = 0; k < param.size(); ++k) {
		    nfactor += paramset.count(param[k]);
		}
		if (nfactor > 1)
		    return false; //quadratic or higher term
		break;
	    default:
		return false;
	    }
	}
	else if (isMixture(dtrm_nodes[j])) {
	    // Check that indices do not depend on snode
	    if (!dtrm_nodes[j]->isLinear(paramset, false)) {
		return false;
	    }
	}
	else {
	    return false;
	}
    }
    */


    return true; //We made it!
}



void ConjugateGamma::update(RNG *rng)
{
    vector<StochasticNode const*> const &stoch_children = stochasticChildren();
    unsigned int nchildren = stoch_children.size();

    double r; // shape
    double mu; // 1/scale
    double xold = *node()->value(chain());

    //Prior
    vector<Node *> const &param = node()->parents();
    switch(_target_dist) {
    case GAMMA:
	r = *param[0]->value(chain());
	mu = *param[1]->value(chain());
	break;
    case EXP:
	r = 1;
	mu = *param[0]->value(chain());
	break;
    case CHISQ:
	r = *param[0]->value(chain())/2;
	mu = 1/2;
	break;
    default:
	throw logic_error("invalid distribution in ConjugateGamma sampler");
    }

    // likelihood 
    bool empty = deterministicChildren().empty();
    bool temp_coef = false;
    if (!empty) {
	if (_coef == 0) {
	    temp_coef = true;
	    _coef = new double[nchildren];
	}
	calCoef();
    }
    for (unsigned int i = 0; i < nchildren; ++i) {

	double coef_i = empty ? 1 : _coef[i];
	if (coef_i > 0) {

	    StochasticNode const *schild = stoch_children[i];
	    vector<Node*> const &cparam = schild->parents();
	    double Y = *schild->value(chain());
	    double ymean; //normal mean
	    switch(_child_dist[i]) {
	    case GAMMA:
		r += *cparam[0]->value(chain());
		mu += coef_i * Y;
		break;
	    case EXP:
		r += 1;
		mu += coef_i * Y;
		break;
	    case NORM:
		r += 0.5;
		ymean = *cparam[0]->value(chain());
		mu += coef_i * (Y - ymean) * (Y - ymean) / 2;
		break;
	    case POIS:
		r += Y;
		mu += coef_i;
		break;
	    case DEXP:
		r += 1;
		ymean = *cparam[0]->value(chain());
		mu += coef_i * fabs(Y - ymean);
		break;
	    case WEIB:
		r += 1; 
		mu += coef_i * pow(Y, *cparam[0]->value(chain()));
		break;
	    default:
		throw logic_error("Invalid distribution in Conjugate Gamma sampler");
	    }
	}
    }
    if (temp_coef) {
	delete [] _coef;
	_coef = 0;
    }

    // Sample from the posterior
    double xnew;
    if (node()->isBounded()) {
	// Use inversion to get random sample
	double lower = 0;
	Node const *lb = node()->lowerBound();
	if (lb) {
	    lower = max(lower, *lb->value(chain()));
	}
	Node const *ub = node()->upperBound();
	double plower = lb ? pgamma(lower,               r, 1/mu, 1, 0) : 0;
	double pupper = ub ? pgamma(*ub->value(chain()), r, 1/mu, 1, 0) : 1;
	double p = runif(plower, pupper, rng);
	xnew = qgamma(p, r, 1/mu, 1, 0);    
    }
    else {
	xnew = rgamma(r, 1/mu, rng);
    }
    setValue(&xnew, 1);  
}
