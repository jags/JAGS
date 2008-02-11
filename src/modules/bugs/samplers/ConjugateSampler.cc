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
using std::string;

static map<string, ConjugateDist> _dist_table;
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

  initialized = true;
}

/* FIXME: _child_dist should be const */

ConjugateDist getDist(StochasticNode const *snode)
{
    if(!initialized)
	init_tables();
  
    string const &name = snode->distribution()->name();
    map<string, ConjugateDist>::iterator p(_dist_table.find(name));

    if (p == _dist_table.end())
	return OTHERDIST;
    else
	return p->second;
}

ConjugateSampler::ConjugateSampler(StochasticNode *node, Graph const &graph, 
				   ConjugateMethod *method)
    : Sampler(vector<StochasticNode*>(1,node), graph),
      _snode(node),
      _method(method),
      _target_dist(getDist(node))
{
    vector<StochasticNode const*> const &children = stochasticChildren();
    for (unsigned int i = 0; i < children.size(); ++i) {
	_child_dist.push_back(getDist(children[i]));
    }
    
    method->initialize(this, graph);
}

ConjugateSampler::~ConjugateSampler()
{
    delete _method;
}

void ConjugateSampler::update(vector<RNG*> const &rngs)
{
    unsigned int nchain = _snode->nchain();
    for (unsigned int ch = 0; ch < nchain; ++ch) {
	_method->update(this, ch, rngs[ch]);
    }
}

bool ConjugateSampler::isAdaptive() const
{
    return false;
}

bool ConjugateSampler::adaptOff()
{
    return true;
}

StochasticNode *ConjugateSampler::node() const
{
    return _snode;
}

vector<ConjugateDist> const &ConjugateSampler::childDist() const
{
    return _child_dist;
}

ConjugateDist ConjugateSampler::targetDist() const
{
    return _target_dist;
}

string ConjugateSampler::name() const
{
    return _method->name();
}
