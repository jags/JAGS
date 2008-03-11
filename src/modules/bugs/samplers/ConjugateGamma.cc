#include <config.h>

#include "ConjugateGamma.h"
#include "ConjugateSampler.h"

#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <graph/MixtureNode.h>
#include <graph/NodeError.h>
#include <sarray/SArray.h>
#include <sampler/Linear.h>

#include <set>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <algorithm>

#include <JRmath.h>

using std::vector;
using std::set;
using std::sqrt;
using std::invalid_argument;
using std::logic_error;
using std::max;
using std::sort;
using std::string;

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

static void calCoef(double *coef, ConjugateSampler *sampler, unsigned int chain)
{   
    const double xold = sampler->node()->value(chain)[0];
    vector<StochasticNode const*> const &stoch_children =
        sampler->stochasticChildren();
    unsigned long nchildren = stoch_children.size();
    vector<ConjugateDist> const &child_dist = sampler->childDist();
    
    for (unsigned int i = 0; i < nchildren; ++i) {
        coef[i] = -getScale(stoch_children[i], child_dist[i], chain);
    }
    double val = xold + 1;
    sampler->setValue(&val, 1, chain);
    for (unsigned int i = 0; i < nchildren; ++i) {
        coef[i] += getScale(stoch_children[i], child_dist[i], chain);
    }
    sampler->setValue(&xold, 1, chain);
}

ConjugateGamma::ConjugateGamma() 
    : _coef(0) 
{
}

void ConjugateGamma::initialize(ConjugateSampler *sampler, Graph const &graph)
{
    // Check for constant scale transformation
    
    if (sampler->deterministicChildren().empty())
	return; //trivial case: all coefficients are 1

    if (!checkScale(sampler->node(), graph, true))
	return; //Not a fixed scale transformation
	
    // Allocate _coef and do one-time calculation of coefficients 
    _coef = new double[sampler->stochasticChildren().size()];
    calCoef(_coef, sampler, 0);
}

ConjugateGamma::~ConjugateGamma()
{
    delete [] _coef;
}


bool ConjugateGamma::canSample(StochasticNode *snode, Graph const &graph)
{
    switch (getDist(snode)) {
    case GAMMA: case EXP: case CHISQ:
	//The exponential and chisquare distributions are both special
	//cases of the gamma distribution and are handled by the conjugate
	//gamma sampler.
	break;
    default:
	return false;
    }

    vector<StochasticNode const*> stoch_nodes;
    vector<Node*> dtrm_nodes;
    Sampler::classifyChildren(vector<StochasticNode*>(1,snode), 
		     graph, stoch_nodes, dtrm_nodes);
    // Create a set of nodes containing snode and its deterministic
    // descendants for the checks below.
    set<Node const *> paramset;
    paramset.insert(snode);
    for (unsigned int j = 0; j < dtrm_nodes.size(); ++j) {
	paramset.insert(dtrm_nodes[j]);
    }

    // Check stochastic children
    for (unsigned int i = 0; i < stoch_nodes.size(); ++i) {
	vector<Node const*> const &param = stoch_nodes[i]->parents();
	if (isBounded(stoch_nodes[i])) {
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
    if (!checkScale(snode, graph, false)) {
	return false;
    }
    return true; //We made it!
}


void ConjugateGamma::update(ConjugateSampler *sampler,
			    unsigned int chain, RNG *rng) const
{
    vector<StochasticNode const*> const &stoch_children = 
	sampler->stochasticChildren();
    unsigned int nchildren = stoch_children.size();

    double r; // shape
    double mu; // 1/scale

    //Prior
    vector<Node const *> const &param = sampler->node()->parents();
    switch(sampler->targetDist()) {
    case GAMMA:
	r = *param[0]->value(chain);
	mu = *param[1]->value(chain);
	break;
    case EXP:
	r = 1;
	mu = *param[0]->value(chain);
	break;
    case CHISQ:
	r = *param[0]->value(chain)/2;
	mu = 1/2;
	break;
    default:
	throw logic_error("invalid distribution in ConjugateGamma sampler");
    }

    // likelihood 
    double *coef = 0;
    bool empty = sampler->deterministicChildren().empty();
    bool temp_coef = false;
    if (!empty && _coef == 0) {
	    temp_coef = true;
	    coef = new double[nchildren];
	    calCoef(coef, sampler, chain);
    }
    else {
	coef = _coef;
    }

    vector<ConjugateDist> const &child_dist = sampler->childDist();
    for (unsigned int i = 0; i < nchildren; ++i) {

	double coef_i = empty ? 1 : coef[i];
	if (coef_i > 0) {

	    StochasticNode const *schild = stoch_children[i];
	    vector<Node const*> const &cparam = schild->parents();
	    double Y = *schild->value(chain);
	    double ymean; //normal mean
	    switch(child_dist[i]) {
	    case GAMMA:
		r += *cparam[0]->value(chain);
		mu += coef_i * Y ;
		break;
	    case EXP:
		r += 1;
		mu += coef_i * Y;
		break;
	    case NORM:
		r += 0.5;
		ymean = *cparam[0]->value(chain);
		mu += coef_i * (Y - ymean) * (Y - ymean) / 2;
		break;
	    case POIS:
		r += Y;
		mu += coef_i;
		break;
	    case DEXP:
		r += 1;
		ymean = *cparam[0]->value(chain);
		mu += coef_i * fabs(Y - ymean);
		break;
	    case WEIB:
		r += 1; 
		mu += coef_i * pow(Y, *cparam[0]->value(chain));
		break;
	    default:
		throw logic_error("Invalid distribution in Conjugate Gamma sampler");
	    }
	}
    }
    if (temp_coef) {
	delete [] coef;
    }

    // Sample from the posterior
    double xnew;
    if (isBounded(sampler->node())) {
	// Use inversion to get random sample
	double lower = 0;
	Node const *lb = sampler->node()->lowerBound();
	if (lb) {
	    lower = max(lower, *lb->value(chain));
	}
	Node const *ub = sampler->node()->upperBound();
	double plower = lb ? pgamma(lower,             r, 1/mu, 1, 0) : 0;
	double pupper = ub ? pgamma(*ub->value(chain), r, 1/mu, 1, 0) : 1;
	double p = runif(plower, pupper, rng);
	xnew = qgamma(p, r, 1/mu, 1, 0);    
    }
    else {
	xnew = rgamma(r, 1/mu, rng);
    }
    sampler->setValue(&xnew, 1, chain);  
}

string ConjugateGamma::name() const
{
    return "ConjugateGamma";
}
