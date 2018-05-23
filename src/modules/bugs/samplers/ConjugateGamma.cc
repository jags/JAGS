#include <config.h>
#include "ConjugateGamma.h"

#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <graph/MixtureNode.h>
#include <graph/NodeError.h>
#include <graph/Graph.h>
#include <sarray/SArray.h>
#include <sampler/Linear.h>
#include <sampler/SingletonGraphView.h>
#include <module/ModuleError.h>

#include <set>
#include <vector>
#include <cmath>
#include <algorithm>

#include <JRmath.h>

using std::vector;
using std::set;
using std::sqrt;
using std::max;
using std::sort;
using std::string;

namespace jags {
namespace bugs {

static double 
getScale(StochasticNode const *snode, ConjugateDist d, unsigned int chain)
{
    //Get scale parameter of snode
    switch(d) {
    case GAMMA: case NORM: case DEXP: case WEIB: case LNORM:
	return *snode->parents()[1]->value(chain);
	break;
    case EXP: case POIS:
	return *snode->parents()[0]->value(chain);
	break;
    default:
	throwNodeError(snode, 
		       "Can't get scale parameter: invalid distribution");
	return 0; //-Wall
	
    } 
}


static void calCoef(double *coef, SingletonGraphView const *gv,
		    vector<ConjugateDist> const &child_dist, unsigned int chain)
{   
    const double xold = gv->node()->value(chain)[0];
    vector<StochasticNode *> const &stoch_children =
        gv->stochasticChildren();
    unsigned long nchildren = stoch_children.size();
    
    for (unsigned long i = 0; i < nchildren; ++i) {
        coef[i] = -getScale(stoch_children[i], child_dist[i], chain);
    }
    double val = xold + 1;
    gv->setValue(&val, 1, chain);
    for (unsigned long i = 0; i < nchildren; ++i) {
        coef[i] += getScale(stoch_children[i], child_dist[i], chain);
    }
    gv->setValue(&xold, 1, chain);
}


ConjugateGamma::ConjugateGamma(SingletonGraphView const *gv)
    : ConjugateMethod(gv), _coef(0)
{
    if(!gv->deterministicChildren().empty() && checkScale(gv, true)) 
    {
	//One-off calculation of fixed scale transformation
	_coef = new double[gv->stochasticChildren().size()];
	calCoef(_coef, gv, _child_dist, 0);
    }

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
	//gamma method.
	break;
    default:
	return false;
    }

    SingletonGraphView gv(snode, graph);

    // Check stochastic children
    vector<StochasticNode *> const &stoch_nodes = 
	gv.stochasticChildren();
    for (unsigned int i = 0; i < stoch_nodes.size(); ++i) {
	if (isBounded(stoch_nodes[i])) {
	    return false; //Bounded
	}
	switch(getDist(stoch_nodes[i])) {
	case EXP: case POIS:
	    break;
	case GAMMA: case NORM: case DEXP: case WEIB: case LNORM:
	    if (gv.isDependent(stoch_nodes[i]->parents()[0])) {
		return false; //non-scale parameter depends on snode
	    }
	    break;
	default:
	    return false;
	}
    }

    // Check deterministic descendants are scale transformations 
    if (!checkScale(&gv, false)) {
	return false;
    }
    return true; //We made it!
}


void ConjugateGamma::update(unsigned int chain, RNG *rng) const
{
    vector<StochasticNode *> const &stoch_children = 
	_gv->stochasticChildren();
    unsigned long nchildren = stoch_children.size();

    //Need to initialize these for -Wall
    double r=0; // shape
    double mu=0; // 1/scale

    //Prior
    vector<Node const *> const &param = _gv->node()->parents();
    switch(_target_dist) {
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
	mu = 0.5; //FIXME: This was wrong in 3.x.y. Needs testing
	break;
    default:
	throwLogicError("invalid distribution in ConjugateGamma method");
    }

    // likelihood 
    double *coef = 0;
    bool empty = _gv->deterministicChildren().empty();
    bool temp_coef = false;
    if (!empty && _coef == 0) {
	    temp_coef = true;
	    coef = new double[nchildren];
	    calCoef(coef, _gv, _child_dist, chain);
    }
    else {
	coef = _coef;
    }

    for (unsigned long i = 0; i < nchildren; ++i) {

	double coef_i = empty ? 1 : coef[i];
	if (coef_i > 0) {

	    StochasticNode const *schild = stoch_children[i];
	    vector<Node const*> const &cparam = schild->parents();
	    double Y = *schild->value(chain);
	    double m = *cparam[0]->value(chain); //location parameter 
	    switch(_child_dist[i]) {
	    case GAMMA:
		r += m;
		mu += coef_i * Y ;
		break;
	    case EXP:
		r += 1;
		mu += coef_i * Y;
		break;
	    case NORM:
		r += 0.5;
		mu += coef_i * (Y - m) * (Y - m) / 2;
		break;
	    case POIS:
		r += Y;
		mu += coef_i;
		break;
	    case DEXP:
		r += 1;
		mu += coef_i * fabs(Y - m);
		break;
	    case WEIB:
		r += 1; 
		mu += coef_i * pow(Y, m);
		break;
	    case LNORM:
		r+= 0.5;
		mu += coef_i * (log(Y) - m) * (log(Y) - m) / 2;
		break;
	    default:
		throwLogicError("Invalid distribution in Conjugate Gamma method"); 
	    }
	}
    }
    if (temp_coef) {
	delete [] coef;
    }

    // Sample from the posterior
    double xnew;
    if (isBounded(_gv->node())) {
	// Use inversion to get random sample
	double lower = 0;
	Node const *lb = _gv->node()->lowerBound();
	if (lb) {
	    lower = max(lower, *lb->value(chain));
	}
	Node const *ub = _gv->node()->upperBound();
	double plower = lb ? pgamma(lower,             r, 1/mu, 1, 0) : 0;
	double pupper = ub ? pgamma(*ub->value(chain), r, 1/mu, 1, 0) : 1;
	double p = runif(plower, pupper, rng);
	xnew = qgamma(p, r, 1/mu, 1, 0);    
    }
    else {
	xnew = rgamma(r, 1/mu, rng);
    }
    _gv->setValue(&xnew, 1, chain);
}

}}
