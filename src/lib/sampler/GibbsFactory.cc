#include <sampler/GibbsFactory.h>
#include <sampler/GibbsSampler.h>
#include <graph/StochasticNode.h>
#include <sampler/Sampler.h>

using std::vector;
using std::set;

void GibbsFactory::makeSampler(set<StochasticNode*> &nodes, 
                               Graph const &graph,
                               vector<vector<Sampler*> > &samplers) const
{
  for (set<StochasticNode*>::iterator p(nodes.begin()); p != nodes.end(); )
    {
      if (canSample(*p, graph)) {
        for (unsigned int n = 0; n < samplers.size(); ++n) {
	   samplers[n].push_back(makeGibbsSampler(*p, graph, n));
        }
	nodes.erase(p++); //Avoid invalidating the iterator
      }
      else {
	++p;
      }
    }
}
