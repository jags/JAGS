#include <config.h>

#include "Censored.h"
#include "DInterval.h"

#include <sampler/SingletonGraphView.h>
#include <graph/NodeError.h>
#include <graph/StochasticNode.h>
#include <graph/Graph.h>
#include <module/ModuleError.h>

#include <vector>
#include <cmath>

using std::vector;
using std::string;

namespace jags {
namespace bugs {

static double indicator(SingletonGraphView const *gv, unsigned int ch)
{
    return gv->stochasticChildren()[0]->value(ch)[0];
}

static Node const *breaks(SingletonGraphView const *gv)
{
    return gv->stochasticChildren()[0]->parents()[1];
}

Censored::Censored(SingletonGraphView const *gv)
    : ConjugateMethod(gv), 
      _snode(dynamic_cast<ScalarStochasticNode*>(gv->node()))
{
    unsigned long nbreaks = breaks(gv)->length();
    for (unsigned int ch = 0; ch < _snode->nchain(); ++ch) {
	double y = indicator(gv, ch);
	if (y < 0 || y > nbreaks) {
	    throwNodeError(_snode, "Bad interval-censored node");
	}
    }
}

bool Censored::canSample(StochasticNode *snode, Graph const &graph)
{
    // The sampler works on scalar real-valued distributions with a
    // single child: an observed stochastic node with a "dinterval"
    // distribution. 
  
    if (snode->isDiscreteValued())
	return false;

    if (!snode->fullRank())
	return false;
    
    if (dynamic_cast<ScalarStochasticNode*>(snode) == nullptr)
	return false;

    // The sampler relies on the fact that boundable distributions can
    // draw truncated random samples. However, the node itself must not
    // be bounded as we don't want the additional complication of combining
    // a priori and a posteriori bounds on the distribution.

    if (!snode->distribution()->canBound())
	return false;
    if (isBounded(snode))
	return false;
  
    //Check that we have a single stochastic child, which is a direct
    //child of the sampled node which has distribution "dinterval"
    SingletonGraphView gv(snode, graph);
    vector<StochasticNode *> const &schild = gv.stochasticChildren();
    vector<DeterministicNode*> const &dchild = gv.deterministicChildren();

    if(schild.size() != 1)
	return false; //Too many children
    if(!dchild.empty())
	return false; //Not direct child
    if (schild[0]->distribution()->name() != "dinterval")
	return false;
    if (schild[0]->parents()[1] == snode)
	return false; //Breaks depend on snode

    return true;
}	

void Censored::update(unsigned int chain, RNG * rng) const
{
    unsigned long y = static_cast<unsigned long>(indicator(_gv, chain));
    double const *b = breaks(_gv)->value(chain);
    unsigned long ymax = breaks(_gv)->length();

    double const *lower = (y == 0) ? nullptr : b + y - 1;
    double const *upper = (y == ymax) ? nullptr : b + y;
	
    _snode->truncatedSample(rng, chain, lower, upper);
}

}}
