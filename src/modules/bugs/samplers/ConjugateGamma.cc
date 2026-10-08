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
#include <util/nainf.h>

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

//debuggin
#include <stdexcept>
#include "SingletonForwardDiff.h"
#include <iostream>

namespace jags {
namespace bugs {

static Node const * 
getScale(StochasticNode const *snode, ConjugateDist d)
{
    //Get scale parameter of snode
    Node const * scale = nullptr;
    switch(d) {
    case GAMMA: case NORM: case DEXP: case WEIB: case LNORM:
	scale = snode->parents()[1];
	break;
    case EXP: case POIS:
	scale = snode->parents()[0];
	break;
    case BERN: case BETA: case BIN: case CAT: case CHISQ: case DIRCH:
    case LOGIS: case MNORM: case MULTI: case NEGBIN: case PAR: case T:
    case UNIF: case WISH: case OTHERDIST:
	throwNodeError(snode,
		       "Can't get scale parameter: invalid distribution");
    }
    return scale;
}

#ifdef DEBUG
static void calCoefNumeric(vector<double> &coef, SingletonGraphView const *gv,
			   vector<ConjugateDist> const &child_dist, unsigned int chain)
{   
    const double xold = gv->node()->value(chain)[0];
    vector<StochasticNode *> const &stoch_children =
        gv->stochasticChildren();
    unsigned long nchildren = stoch_children.size();
    
    for (unsigned long i = 0; i < nchildren; ++i) {
        coef[i] = -getScale(stoch_children[i], child_dist[i])->value(chain)[0];
    }
    double val = xold + 1;
    gv->setValue(&val, 1, chain);
    for (unsigned long i = 0; i < nchildren; ++i) {
        coef[i] += getScale(stoch_children[i], child_dist[i])->value(chain)[0];
    }
    gv->setValue(&xold, 1, chain);
}
#endif

ConjugateGamma::ConjugateGamma(SingletonGraphView const *gv)
    : ConjugateMethod(gv), _fd(nullptr), _coef(nullptr)
{
    vector<StochasticNode *> const &stoch_children = gv->stochasticChildren();
    unsigned long nchildren = stoch_children.size();
    vector<Node const *> dependent_nodes(nchildren);
    
    for (unsigned long i = 0; i < nchildren; ++i) {
	dependent_nodes.push_back(getScale(stoch_children[i], _child_dist[i]));
    }
    _fd = new SingletonForwardDiff(gv, dependent_nodes);
    
    if(!gv->deterministicChildren().empty() && checkScale(gv, true)) 
    {
	//One-off calculation of fixed scale transformation
	_coef = new double[gv->stochasticChildren().size()];
	_fd->calCoef(_coef, 0);
	delete _fd;
	_fd = nullptr;

#ifdef DEBUG
	vector<double> betas_num(gv->length());
	calCoefNumeric(betas_num, gv, _child_dist, 0);

	for (unsigned int i = 0; i < betas_num.size(); ++i) {
	    if (fabs(_coef[i] - betas_num[i]) > 1.0E-6) {
		throw std::logic_error("ConjugateGamma: Error in scale calculation");
	    }
	}
#endif
    }

}

ConjugateGamma::~ConjugateGamma()
{
    delete [] _coef;
    delete _fd;
}


bool ConjugateGamma::canSample(StochasticNode *snode, Graph const &graph)
{
    switch (getDist(snode)) {
    case GAMMA: case EXP: case CHISQ:
	//The exponential and chisquare distributions are both special
	//cases of the gamma distribution and are handled by the conjugate
	//gamma method.
	break;
    case BERN: case BETA: case BIN: case CAT: case DEXP: case DIRCH:
    case LNORM: case LOGIS: case MNORM: case MULTI: case NEGBIN: case NORM:
    case PAR: case POIS: case T: case UNIF: case WEIB: case WISH:
    case OTHERDIST:
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
	case BERN: case BETA: case BIN: case CAT: case CHISQ: case DIRCH:
	case LOGIS: case MNORM: case MULTI: case NEGBIN: case PAR: case T:
	case UNIF: case WISH: case OTHERDIST:	    
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
    case BERN: case BETA: case BIN: case CAT: case DEXP: case DIRCH:
    case LNORM: case LOGIS: case MNORM: case MULTI: case NEGBIN: case NORM:
    case PAR: case POIS: case T: case UNIF: case WEIB: case WISH:
    case OTHERDIST:	
	throwLogicError("invalid distribution in ConjugateGamma method");
    }

    // likelihood 
    double *coef = nullptr;
    bool empty = _gv->deterministicChildren().empty();
    bool temp_coef = false;
    if (!empty && _coef == nullptr) {
	temp_coef = true;
	coef = new double[nchildren];
	_fd->calCoef(coef, chain);
	    
#ifdef DEBUG
	vector<double> coef_num(nchildren);
	calCoefNumeric(coef_num, _gv, _child_dist, chain);
	for (unsigned int i = 0; i < coef_num.size(); ++i) {
	    if (fabs(coef[i] - coef_num[i]) > 1.0E-6) {
		throwLogicError("ConjugateGamma: Error in scale calculation");
	    }
	}
#endif
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
	    case BERN: case BETA: case BIN: case CAT: case CHISQ: case DIRCH:
	    case LOGIS: case MNORM: case MULTI: case NEGBIN: case PAR: case T:
	    case UNIF: case WISH: case OTHERDIST:	
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
