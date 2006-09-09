#include <config.h>
#include <graph/StochasticNode.h>
#include <graph/LogicalNode.h>
#include <distribution/Distribution.h>
#include <function/Function.h>

#include <map>
#include <string>

#include "ConjugateSampler.h"

using std::map;
using std::string;
using std::vector;

static map<const string, ConjugateDist> _dist_table;
//static map<const string, ConjugateOp> _op_table;
static bool initialized = false;

static void init_tables ()
{
  _dist_table["dbern"] = BERN;
  _dist_table["dbeta"] = BETA;
  _dist_table["dbin"] = BIN;
  _dist_table["dcat"] = CAT;
  _dist_table["dchisq"] = CHISQ;
  _dist_table["ddexp"] = DEXP;
  _dist_table["ddirch"] = DIRCH;
  _dist_table["dexp"] = EXP;
  _dist_table["dgamma"] = GAMMA;
  _dist_table["dlnorm"] = LNORM;
  _dist_table["dlogis"] = LOGIS;
  _dist_table["dmnorm"] = MNORM;
  _dist_table["dmulti"] = MULTI;
  _dist_table["dnegbin"] = NEGBIN;
  _dist_table["dnorm"] = NORM;
  _dist_table["dpar"] = PAR;
  _dist_table["dpois"] = POIS;
  _dist_table["dt"] = T;
  _dist_table["dunif"] = UNIF;
  _dist_table["dweib"] = WEIB;
  _dist_table["dwish"] = WISH;

  /*
  _op_table["*"] = MULTIPLY;
  _op_table["/"] = DIVIDE;
  _op_table["+"] = ADD;
  _op_table["-"] = SUBTRACT;
  _op_table["NEG"] = NEG;
  _op_table["inprod"] = INPROD;
  */

  initialized = true;
}


ConjugateDist getDist(StochasticNode const *snode)
{
  if(!initialized)
    init_tables();
  
  string const &name = snode->distribution()->name();
  map<const string, ConjugateDist>::iterator p(_dist_table.find(name));

  if (p == _dist_table.end())
    return OTHERDIST;
  else
    return p->second;
}

/*
ConjugateOp getOp(LogicalNode const *lnode)
{
  if(!initialized)
    init_tables();
  
  string const &name = lnode->function()->name();
  map<const string, ConjugateOp>::iterator p(_op_table.find(name));
  if (p == _op_table.end())
    return OTHEROP;
  else
    return p->second;
}
*/

ConjugateSampler::ConjugateSampler(StochasticNode *node, Graph const &graph, unsigned int chain)
  : GibbsSampler(node, graph, chain), _target_dist(getDist(node))
{
  vector<StochasticNode const*> const &children = stochasticChildren();
  for (unsigned int i = 0; i < children.size(); ++i) {
    _child_dist.push_back(getDist(children[i]));
  }
}

void ConjugateSampler::burninOff()
{
}
