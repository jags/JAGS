#include <config.h>

#include "Censored.h"
#include "DInterval.h"
//#include "RNG.h"

#include <graph/NodeError.h>
#include <graph/StochasticNode.h>
#include <distribution/DistReal.h>

#include <stdexcept>
#include <vector>
#include <cmath>

#include <Rmath.h>

using std::vector;
using std::invalid_argument;

Censored::Censored(StochasticNode *snode, Graph const &graph, unsigned int chain)
  : GibbsSampler(snode, graph, chain)
{
  if (!canSample(snode, graph)) {
    throw invalid_argument("Can't construct Censored sampler");
  }
}

Censored::~Censored()
{
}

bool Censored::canSample(StochasticNode *snode, Graph const &graph)
{
  /* 
     FIXME: Doesn't currently work without downcasting to incomplete
     type, so for the moment just disable it.
  */
  return false;

  /* The sampler works on scalar real-valued distributions with a
     single child: an observed stochastic node with a "dinterval"
     distribution
  */
  
  if(snode->distribution()->isDiscreteValued() || snode->length() != 1)
    return false;
  
  //FIXME: In principle, we could allow bounded nodes.
  if(snode->isBounded())
    return false;
  
  vector<StochasticNode const*> stoch_nodes;
  vector<DeterministicNode*> dtrm_nodes;
  classifyChildren(vector<StochasticNode*>(1,snode), 
		   graph, stoch_nodes, dtrm_nodes);

  if (!dtrm_nodes.empty())
    return false;

  if (stoch_nodes.size() == 1) {
    return (stoch_nodes[0]->distribution()->name() == "dinterval");
  }
  else {
    return false;
  }
}

void Censored::update(double (*unif_rand)(void))
{
  /*
  StochasticNode const *child = stochasticChildren()[0];
  double const *breaks = child->parameters()[1]->value();
  unsigned long nbreak = child->parameters()[1]->length();
  long y = static_cast<long>(child->data.value()[0]);  

  if (y < 0 || y > nbreak)
    throw NodeError(node(), "Bad interval-censored node");

  Distribution dist = node()->distribution();
  vector<SArray const *> const &parameters = node()->parameters();

  double x;
  if (y==0) {
    //Left-censored: x <= breaks[0]
    double logpx = dreal->p(breaks[0], parameters, true, true) - rexp(1);
    x = dreal->q(logpx, parameters, true, true);
  }
  else if (y==nbreak) {
    //Right-censored: x > breaks[nbreak-1]
    double logpx = dreal->p(breaks[nbreak-1], parameters, false, true) 
      - rexp(1);
    x = dreal->q(logpx, parameters, false, true);
  }
  else {
    //Interval-censored: breaks[y-1] < x <= breaks[y]
    double px = runif(breaks[y-1], breaks[y]);
    x = dreal->q(px, parameters, true, false);
  }

  setValue(&x,1);
  */
}

void Censored::burninOff()
{}
