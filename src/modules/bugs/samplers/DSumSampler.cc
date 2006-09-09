#include <config.h>
#include <distribution/DistDiscrete.h>
#include <graph/StochasticNode.h>
#include <graph/Graph.h>
#include "DSumSampler.h"

#include <cfloat>
#include <climits>
#include <stdexcept>
#include <set>
#include <algorithm>
#include <cmath>

#include <Rmath.h>

using std::set;
using std::vector;
using std::invalid_argument;
using std::max;
using std::min;
using std::exp;

DSumSampler::DSumSampler(vector<StochasticNode *> const &nodes,
			 Graph const &graph, unsigned int chain)
    : Slicer(nodes, graph, chain, 2, 10)
{
    if (!canSample(nodes, graph)) {
	throw invalid_argument("Can't construct DSumSampler");
    }

    Node const *dsum = 0;
    set<Node*> const &children = nodes[0]->children();
    for (set<Node*>::const_iterator p = children.begin(); p != children.end();
	 ++p) 
	{
	  if (asStochastic(*p) && (*p)->isObserved() &&
	      asStochastic(*p)->distribution()->name() == "dsum")
	    {
	      dsum = *p;
	      break;
	    }
	}
  
    _sum = static_cast<long>(*dsum->value(chain));
    _x = *nodes[0]->value(chain);
    //Make sure values are consistent at start
    double x2 = _sum - static_cast<long>(_x);
    nodes[1]->setValue(&x2,1,chain);
}

DSumSampler::~DSumSampler()
{
}

bool DSumSampler::canSample(vector<StochasticNode *> const &nodes,
			    Graph const &graph)
{
  if (nodes.size() != 2)
    return false;

  if (!(graph.contains(nodes[0]) && graph.contains(nodes[1])))
    return false;

  Node const *dsum[2] = {0,0};
  for (unsigned int i = 0; i < 2; ++i) {
    // Nodes must be scalar ...
    if (nodes[i]->length() != 1)
      return false;
    
    // stochastic ...
    if (!asStochastic(nodes[i]))
      return false;

    // discrete-valued ...
    if (!nodes[i]->isDiscreteValued())
      return false;

    // Nodes must be parents of an observed dsum node
    set<Node*> const &children = nodes[i]->children();
    for (set<Node*>::const_iterator p = children.begin(); 
	 p != children.end(); ++p) 
      {
	if (asStochastic(*p) && (*p)->isObserved() &&
	    asStochastic(*p)->distribution()->name() == "dsum")
	  {
	    dsum[i] = *p;
	    break;
	  }
      }
    if (dsum[i] == 0) {
      return false;
    }
  }
  // Must be the same dsum node!
  if (dsum[0] != dsum[1])
    return false;

  return true;
}

/*
void DSumSampler::update()
{
  // Calculate lower and upper limits of first variable
  DistDiscrete const *idist0 = asDiscrete(nodes()[0]->distribution());
  DistDiscrete const *idist1 = asDiscrete(nodes()[1]->distribution());
  long ll = max(idist0->lowerLimit(nodes()[0]->parameters()),
		_sum - idist1->upperLimit(nodes()[1]->parameters()));
  long ul = min(idist0->upperLimit(nodes()[0]->parameters()),
		_sum - idist1->lowerLimit(nodes()[1]->parameters()));

  long N = ul - ll + 1;
  double *p = new double[N];
  
  //Setup the vectors used by Sampler::setValue
  double v[2];
  vector<double const*> value;
  value.push_back(&v[0]);
  value.push_back(&v[1]);
  vector<long> length = vector<long>(2,1);

  //Calculate likelihood for each possible value
  for (int i = 0; i < N; ++i) {
    v[0] = ll + i;
    v[1] = _sum - ll - i;
    setValue(value, length);
    p[i] = exp(logFullConditional());
  }

  //Normalize
  double sump = 0;
  for (int i = 0; i < N; ++i) {
    sump += p[i];
  }

  //Sample
  double urand = runif(0.0, sump);
  long i;
  sump = 0;
  for (i = 0; i < N - 1; ++i) {
    sump += p[i];
    if (sump > urand) {
      break;
    }
  }

  //Set sampled value
  v[0] = ll + i;
  v[1] = _sum - ll - i;
  setValue(value, length);

  delete [] p;
}
*/

void DSumSampler::setValue(double x)
{
  double v[2];
  vector<double const*> value;
  value.push_back(&v[0]);
  value.push_back(&v[1]);
  vector<unsigned long> length = vector<unsigned long>(2,1);  

  _x = x;
  v[0] = static_cast<long>(x);
  v[1] = _sum - v[0];
 
  Sampler::setValue(value, length);
}

double DSumSampler::value()
{
  return _x;
}

double DSumSampler::lowerLimit()
{
  vector<StochasticNode *> const &n = nodes();
  double l0 = n[0]->distribution()->lowerSupport(0, n[0]->parameters(chain()));
  double u1 = n[1]->distribution()->upperSupport(0, n[1]->parameters(chain()));
  return max(l0, _sum - u1);
}

double DSumSampler::upperLimit()
{
  vector<StochasticNode *> const &n = nodes();
  double u1 = n[0]->distribution()->upperSupport(0, n[0]->parameters(chain()));
  double l2 = n[1]->distribution()->lowerSupport(0, n[1]->parameters(chain()));
  return min(u1, _sum - l2);
}

void DSumSampler::update(RNG *rng)
{
  updateDouble(rng);
}
