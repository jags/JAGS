#include <config.h>

#include "TruncatedGamma.h"

#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <graph/MixtureNode.h>
#include <graph/NodeError.h>
#include <sarray/SArray.h>
#include <sampler/Linear.h>
#include <sampler/SingletonGraphView.h>
#include <util/nainf.h>
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

#define TG_EPS 0.1

/* 
   The machinery of the TruncatedGamma method is very similar to the
   Conjugate gamma method, except that Poisson distributions are not
   permitted.  
*/

namespace jags {
namespace bugs {

static Node const * getParent(StochasticNode const *snode)
{
    //Choose appropriate parent for stochastic child

    switch(getDist(snode)) {
    case GAMMA: case NORM: case DEXP: case WEIB: case LNORM:
	return snode->parents()[1];
	break;
    case EXP: case POIS:
	return snode->parents()[0];
	break;
    default:
        throwLogicError("Invalid distribution in TruncatedGamma sampler");
	return 0; //-Wall
    } 
}

static double calExponent(SingletonGraphView const *gv, unsigned long i, 
			  unsigned int chain)
{   
    /*        
       All stochastic children depend on the sampled node "x" via a
       parent that takes the value "a * x^c". This function calculates
       the value of the exponent c for the ith stochastic child
    */

    if (gv->deterministicChildren().empty()) {
	return 1;
    }

    //Save original value
    const double xold = *gv->node()->value(chain);

    double x0 = xold;
    if (x0 <= 0) {
	x0 = 1;
	gv->setValue(&x0, 1, chain);
    }
    Node const *par = getParent(gv->stochasticChildren()[i]);
    double y0 = *par->value(chain);
    double x1 = 2 * x0;
    gv->setValue(&x1, 1, chain);
    double y1 = *par->value(chain);

    gv->setValue(&xold, 1, chain); //Restore original value
    
    if (y1 <= 0 || y0 <= 0)
	throwLogicError("Invalid scale function in TruncatedGamma sampler");

    return (log(y1) - log(y0))/log(2.0);
}

TruncatedGamma::TruncatedGamma(SingletonGraphView const *gv)
    : ConjugateMethod(gv), _exponent(calExponent(gv, 0, 0))
{
}

TruncatedGamma::~TruncatedGamma()
{
}

bool TruncatedGamma::canSample(StochasticNode *snode, Graph const &graph)
{
    // Target node must have a uniform distribution with a fixed,
    // non-negative lower bound.
    if (getDist(snode) != UNIF)
	return false;
    if (!snode->parents()[0]->isFixed())
	return false;
    if (snode->parents()[0]->value(0)[0] < 0)
	return false;

    SingletonGraphView gv(snode, graph);
    vector<StochasticNode *> const &stoch_nodes = gv.stochasticChildren();

    //Check that stochastic children are valid
    for (unsigned int i = 0; i < stoch_nodes.size(); ++i) {
	//Check that children are unbounded
	if (isBounded(stoch_nodes[i])) {
	    return false; 
	}
	//Check that are no unwanted paths to the stochastic children 
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

    // Deterministic children must be fixed power functions
    if (!checkPower(&gv, true))	return false;

    // Calculate the coefficient of the power function
    double exp0 = calExponent(&gv, 0, 0);
    for (unsigned long i = 0; i < stoch_nodes.size(); ++i) {
	double expi = calExponent(&gv, i, 0);
	//Exclude exponents near 0 for numerical stability
	if (fabs(expi) < TG_EPS) return false;
	// Exponent must be the same for all stochastic children
	// We allow for some numerical error
	if (fabs(expi - exp0) > 1E-16) return false;
    }

    //A negative prior value of shape is required if we can guarantee
    //that it is positive a posteriori.
    double shape = 1/exp0; 
    for (unsigned int i = 0; i < stoch_nodes.size(); ++i) {
	switch(getDist(stoch_nodes[i])) {
	case GAMMA:
	    if (stoch_nodes[i]->parents()[0]->isFixed()) {
		shape += *stoch_nodes[0]->parents()[0]->value(0);
	    }
	    break;
	case POIS:
	    if (stoch_nodes[i]->isFixed()) {
		shape += *stoch_nodes[i]->value(0);
	    }
	    break;
	case EXP: case DEXP: case WEIB:
	    shape += 1;
	    break;
	case NORM: case LNORM:
	    shape += 0.5;
	    break;
	default:
	    return false;
	}
    }

    // We could allow any shape > 0, but we add a fudge factor to 
    // guarantee stability
    return (shape >= TG_EPS);
}


void TruncatedGamma::update(unsigned int chain, RNG *rng) const
{
    //If "x" is the value of the sampled node, the parameter "z =
    //x^_exponent" has a truncated gamma distribution.  We work in
    //terms of "z" and then translate back to "x"

    //Prior parameters
    double r = 1/_exponent; // shape
    double mu = 0; // 1/scale

    //Save current value
    StochasticNode *snode = _gv->node();
    const double xold = snode->value(chain)[0];
    const double zold = pow(xold, _exponent);

    //Add likelihood contributions to shape and scale parameters
    vector<StochasticNode *> const &schildren = _gv->stochasticChildren();
    for (unsigned int i = 0; i < schildren.size(); ++i) {

	double Y = *schildren[i]->value(chain);
	double m = *schildren[i]->parents()[0]->value(chain); 
	double c = *getParent(schildren[i])->value(chain) / zold;
	if (c > 0) {
	    switch(_child_dist[i]) {
	    case GAMMA:
		r += m;
		mu += c * Y ;
		break;
	    case EXP:
		r += 1;
		mu += c * Y;
		break;
	    case NORM:
		r += 0.5;
		mu += c * (Y - m) * (Y - m) / 2;
		break;
	    case POIS:
		r += Y;
		mu += c;
		break;
	    case DEXP:
		r += 1;
		mu += c * fabs(Y - m);
		break;
	    case WEIB:
		r += 1; 
		mu += c * pow(Y, m);
		break;
	    case LNORM:
		r += 0.5;
		mu += c * (log(Y) - m) * (log(Y) - m) / 2;
		break;
	    default:
		throwLogicError("Invalid distribution in TruncatedGamma");
	    }
	}
    }
    if (mu == 0) {
	
	throwNodeError(snode, "Degenerate posterior in TruncatedGamma sampler");
    }

    //Find the boundaries on the scale of x
    double lx = snode->parents()[0]->value(chain)[0];
    double ux = snode->parents()[1]->value(chain)[0];
    if (xold < lx || xold > ux) {
	throwLogicError("Current value invalid in TruncatedGamma method");
    }


    // Find boundaries on the scale of z
    double ly, uy;
    if (_exponent > 0) {
	ly = (lx > 0) ? exp(_exponent * log(lx)) : 0;
	uy = exp(_exponent * log(ux));
    }
    else {
	// Boundaries are inverted
	ly = exp(_exponent * log(ux));
	uy = (lx > 0) ? exp(_exponent * log(lx)) : JAGS_POSINF;
    }

    double plower = pgamma(ly, r, 1/mu, 1, 0);
    double pupper = pgamma(uy, r, 1/mu, 1, 0);

    // Sample from the posterior
    double ynew;
    if (pupper - plower > 0.5) {
	do {
	    ynew = rgamma(r, 1/mu, rng);
	}
	while (ynew < ly || ynew > uy);
    }
    else {
	double p = runif(plower, pupper, rng);
	ynew = qgamma(p, r, 1/mu, 1, 0);    
    }

    double xnew = exp(log(ynew)/_exponent);
    _gv->setValue(&xnew, 1, chain);  
}

}}
