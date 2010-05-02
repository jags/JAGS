#include <config.h>

#include "TruncatedGamma.h"

#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <graph/MixtureNode.h>
#include <graph/NodeError.h>
#include <sarray/SArray.h>
#include <sampler/Linear.h>
#include <sampler/GraphView.h>
#include <util/nainf.h>

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

/* 
   The machinery of the TruncatedGamma method is very similar to the
   Conjugate gamma method, except that Poisson distributions are not
   permitted.  
*/

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
	return 0;
    } 
}

static double calShape(GraphView const *gv, unsigned int chain)
{   
    //Calculate prior shape parameter

    if (gv->deterministicChildren().empty()) {
	return 1;
    }

    //Save original value
    const double xold = gv->nodes()[0]->value(chain)[0];

    double x0 = xold;
    if (x0 <= 0) {
	x0 = 1;
    }
    Node const *par = getParent(gv->stochasticChildren()[0]);
    double y0 = par->value(chain)[0];
    double x1 = 2 * x0;
    gv->setValue(&x1, 1, 0);
    double y1 = par->value(chain)[0];

    //Restore original value
    gv->setValue(&xold, 1, chain);
    
    if (y1 <= 0 || y0 <= 0)
	throw logic_error("Invalid scale function in TruncatedGamma");

    return log(2.0)/(log(y1) - log(y0));
}

TruncatedGamma::TruncatedGamma(GraphView const *gv)
    : ConjugateMethod(gv), _shape(calShape(gv,0)),
      _parent(getParent(gv->stochasticChildren()[0]))
{
}

TruncatedGamma::~TruncatedGamma()
{
}

bool TruncatedGamma::canSample(StochasticNode *snode, Graph const &graph)
{
    /* Target node must have a uniform distribution with a fixed,
       non-negative lower bound
    */
    if (getDist(snode) != UNIF)
	return false;
    if (!snode->parents()[0]->isObserved())
	return false;
    if (snode->parents()[0]->value(0)[0] < 0)
	return false;

    GraphView gv(vector<StochasticNode*>(1,snode), graph);

    /* Deterministic children must be fixed power functions */
    if (!checkPower(&gv, true)) {
	return false;
    }

    /* Power function must be increasing */
    double shape = calShape(&gv, 0);

    // Check stochastic children
    vector<StochasticNode const*> const &stoch_nodes = gv.stochasticChildren();
    Node const *par0 = getParent(stoch_nodes[0]);
    for (unsigned int i = 0; i < stoch_nodes.size(); ++i) {

	//Check that stochastic children are not bounded
	if (isBounded(stoch_nodes[i])) {
	    return false;
	}
	
	//Check that parent is common to all stochastic children
	Node const *pari = getParent(stoch_nodes[i]);
	if (pari == 0 || pari != par0) {
	    return false;
	}

	ConjugateDist dist = getDist(stoch_nodes[i]);
	vector<Node const*> const &cparam = stoch_nodes[i]->parents();

	//Check that there are no unwanted paths to the stochastic children

	switch(dist) {
	case EXP: case POIS:
	    break;
	case GAMMA: case NORM: case DEXP: case WEIB: case LNORM:
	    if (gv.isDependent(cparam[0])) {
		return false; //non-scale parameter depends on snode
	    }
	    break;
	default:
	    return false;
	}


	//Check that we can guarantee a positive shape parameter a posteriori
	switch(dist) {
	case GAMMA:
	    if (cparam[0]->isObserved()) {
		shape += *cparam[0]->value(0);
	    }
	    break;
	case POIS:
	    if (stoch_nodes[i]->isObserved()) {
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

    /* We could allow any shape > 0, but we add a fudge factor to 
     * guarantee stability
     */
    return (shape >= 0.1);
}


void TruncatedGamma::update(unsigned int chain, RNG *rng) const
{
    //Prior
    double r = _shape; // shape
    double mu = 0; // 1/scale

    // likelihood 
    vector<StochasticNode const*> const &schildren = _gv->stochasticChildren();
    for (unsigned int i = 0; i < schildren.size(); ++i) {

	double Y = *schildren[i]->value(chain);
	double m = *schildren[i]->parents()[0]->value(chain); //location param 
	switch(_child_dist[i]) {
	case GAMMA:
	    r += m;
	    mu += Y ;
	    break;
	case EXP:
	    r += 1;
	    mu += Y;
	    break;
	case NORM:
	    r += 0.5;
	    mu += (Y - m) * (Y - m) / 2;
	    break;
	case POIS:
	    r += Y;
	    mu += 1;
	    break;
	case DEXP:
	    r += 1;
	    mu += fabs(Y - m);
	    break;
	case WEIB:
	    r += 1; 
	    mu += pow(Y, m);
	    break;
	case LNORM:
	    r += 0.5;
	    mu += (log(Y) - m) * (log(Y) - m) / 2;
	    break;
	default:
	    throw logic_error("Invalid distribution in TruncatedGamma");
	}
    }

    //FIXME: This is a hack for badly initialized models. We would not
    //need it if we sampled nodes in the correct order - 2 May 2010
    if (mu == 0)
	return;

    //Save current value
    StochasticNode *snode = _gv->nodes()[0];
    double xold = snode->value(chain)[0];
    double yold = _parent->value(chain)[0];

    //Find the boundaries on the scale of x
    double lx = snode->parents()[0]->value(chain)[0];
    double ux = snode->parents()[1]->value(chain)[0];
    if (xold < lx || xold > ux) {
	throw logic_error("Current value invalid TruncatedGamma method");
    }
    if (lx > ux) {
	throw logic_error("Inconsistent prior in TruncatedGamma method");
    }


    // Find boundaries on the scale of y
    // We need to take care that they may be inverted

    double ly = 0;
    double uy = JAGS_POSINF;

    if (_shape > 0) {

	if (lx > 0) {

	    ly = yold * exp((log(lx) - log(xold))/_shape);
	}
	uy = yold * exp((log(ux) - log(xold))/_shape);

    }
    else {

	if (lx > 0) {
	    uy = yold * exp((log(lx) - log(xold))/_shape);
	}
	ly = yold * exp((log(ux) - log(xold))/_shape);
    }

    double plower = pgamma(ly, r, 1/mu, 1, 0);
    double pupper = pgamma(uy, r, 1/mu, 1, 0);

    // Sample from the posterior
    double ynew;
    if ((pupper - plower) > 0.5) {
	do {
	    ynew = rgamma(r, 1/mu, rng);
	}
	while (ynew < ly || ynew > uy);
    }
    else {
	double p = runif(plower, pupper, rng);
	ynew = qgamma(p, r, 1/mu, 1, 0);    
    }

    double xnew = xold * exp(_shape * (log(ynew) - log(yold)));
    _gv->setValue(&xnew, 1, chain);  
}

string TruncatedGamma::name() const
{
    return "TruncatedGamma";
}
