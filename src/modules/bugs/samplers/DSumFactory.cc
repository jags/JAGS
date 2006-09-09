#include <config.h>
#include <distribution/Distribution.h>
//asStochastic
#include <graph/StochasticNode.h>
//asDiscrete
#include <distribution/DistDiscrete.h>
#include <graph/Graph.h>

#include "DSumFactory.h"
#include "DSumSampler.h"

#include <stdexcept>

using std::set;
using std::vector;
using std::runtime_error;

static bool isDSumNode(Node const *node)
{
  StochasticNode const *snode = asStochastic(node);
  if (snode) {
    return snode->distribution()->name() == "dsum";
  }
  else {
    return false;
  }
}

static bool isDiscreteStochastic(Node const *node)
{
  if (!asStochastic(node))
    return false;
  else if (node->length() != 1)
    return false;
  else
    return node->isDiscreteValued();
}

static bool canSample(vector<Node*> const &parameters, Graph const &graph)
{
  for (vector<Node*>::const_iterator p = parameters.begin();
       p != parameters.end(); ++p)
    {
      Node const *param = *p;
      
      if (!isDiscreteStochastic(param))
	return false;

      /*
	Check that there is only a single child of param within the graph, 
	and that it is a dsum node
      */
      int nchild = 0;
      for (set<Node*>::const_iterator i = param->children().begin();
	   i != param->children().end(); ++i)
	{
	  if (!isDSumNode(*i))
	    return false;

	  ++nchild;
	}
      if (nchild != 1)
	return false;
    }
  return true;
}

void DSumFactory::makeSampler(set<StochasticNode*> &nodes,
			      Graph const &graph,
			      vector<vector<Sampler*> > &samplers) const
{
  set<StochasticNode const*> dsum_nodes;

  for (set<StochasticNode*>::iterator p = nodes.begin(); p != nodes.end(); ++p)
    {
      // Find integer-valued stochastic nodes in the graph ...
      if (isDiscreteStochastic(*p) && graph.contains(*p)) {
	set<Node*> const &children = (*p)->children();
	for (set<Node*>::iterator q = children.begin() ; q != children.end(); 
	     ++q) 
	  {
	    // ... that have a DSum node as a child in the graph
	    if (isDSumNode(*q) && graph.contains(*q)) {
	      dsum_nodes.insert(asStochastic(*q));
	    }
	  }
      }
    }
  
  if (dsum_nodes.empty())
    return; //Nothing to do
  
  for (set<StochasticNode const*>::iterator p = dsum_nodes.begin(); 
       p != dsum_nodes.end(); ++p) 
    {
      vector<Node*> const &parents = (*p)->parents();
      if (canSample(parents, graph)) {
	vector<StochasticNode*> parameters;
	for (vector<Node*>::const_iterator q = parents.begin(); q != parents.end(); 
	     ++q) 
	  {
	    StochasticNode *snode = dynamic_cast<StochasticNode*>(*q);
	    parameters.push_back(snode);
	    nodes.erase(snode);
	  }
        for (unsigned int ch = 0; ch < samplers.size(); ++ch) {
	   samplers[ch].push_back(new DSumSampler(parameters, graph, ch));
        }
      }
      else {
	throw runtime_error("Illegal use of dsum node");
      }
    }
}
