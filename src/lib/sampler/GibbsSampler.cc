#include <config.h>
#include <sampler/GibbsSampler.h>
#include <graph/StochasticNode.h>
#include <graph/DeterministicNode.h>

#include <vector>
#include <set>

using std::vector;
using std::set;

GibbsSampler::GibbsSampler(StochasticNode *node, 
			   Graph const &graph, unsigned int chain)
    : Sampler(vector<StochasticNode *>(1, node), graph, chain)
{
}

void GibbsSampler::setValue(double const *value, unsigned long length)
{
  nodes().front()->setValue(value, length, chain());
  vector<DeterministicNode*> const &dc = deterministicChildren();
  for (vector<DeterministicNode*>::const_iterator i = dc.begin(); 
       i != dc.end(); ++i) 
    {
      (*i)->deterministicSample(chain());
    }
}

StochasticNode *GibbsSampler::node()
{
  return nodes().front();
}
