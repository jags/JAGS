#include <config.h>
#include <sampler/Sampler.h>
#include <graph/StochasticNode.h>
#include <graph/DeterministicNode.h>
#include <graph/Graph.h>
#include <graph/NodeError.h>

#include <stdexcept>
#include <set>
#include <string>
#include <cfloat>
#include <cmath>

using std::vector;
using std::set;
using std::runtime_error;
using std::logic_error;
using std::string;

Sampler::Sampler(vector<StochasticNode *> const &nodes, Graph const &graph,
                 unsigned int chain)
  : _chain(chain), _nodes(nodes)
{
  classifyChildren(nodes, graph, _stoch_children, _determ_children);
}

Sampler::~Sampler()
{
  // virtual destructor
}

vector<StochasticNode *> const &Sampler::nodes() const
{
  return _nodes;
}

static void classifyNode(Node *node, Graph const &sample_graph, 
			 Graph &sgraph, Graph &dgraph)
{
    if (!sample_graph.contains(node))
	return;

    if (node->isVariable() && asStochastic(node)) {
	sgraph.add(node);
    }
    else if (!dgraph.contains(node)) {
	dgraph.add(node);      
	for (set<Node*>::iterator p = node->children().begin();
	     p != node->children().end(); ++p) {
	    classifyNode(*p, sample_graph, sgraph, dgraph);
	}
    }
}

void Sampler::classifyChildren(vector<StochasticNode *> const &nodes,
			       Graph const &graph,
			       vector<StochasticNode const*> &stoch_nodes,
			       vector<DeterministicNode*> &dtrm_nodes)
{
  Graph dgraph, sgraph;

  /* Classify children of each node */
  vector<StochasticNode  *>::const_iterator p = nodes.begin();
  for (; p != nodes.end(); ++p) {
    StochasticNode const *snode = *p;
    if (!graph.contains(snode)) {
      throw logic_error("Sampled node outside of sampling graph");
    }
    for (set<Node*>::const_iterator q = snode->children().begin(); 
	 q != snode->children().end(); ++q) 
      {
	classifyNode(*q, graph, sgraph, dgraph);
      }
  }

  /* Strip nodes to be sampled out of the graph of stochastic
     children. Such nodes would contribute to both the prior
     AND the likelihood, causing incorrect calculation of the
     log full conditional */
  for (p = nodes.begin(); p != nodes.end(); ++p) {
    sgraph.remove(*p);
  }

  vector<Node*> svector;
  sgraph.getNodes(svector);
  for (vector<Node*>::iterator i = svector.begin(); i != svector.end(); 
       ++i) 
    {
      stoch_nodes.push_back(asStochastic(*i));
    }
  
  vector<Node*> dvector;
  dgraph.getSortedNodes(dvector);
  for (vector<Node*>::iterator i = dvector.begin(); i != dvector.end();
       ++i)
    {
      dtrm_nodes.push_back(dynamic_cast<DeterministicNode*>(*i));
    }
}

double Sampler::logFullConditional()
{
  /* We don't know if isinf is a macro or a function declared in the
     standard name space */
  using namespace std;

  double logprior = 0;
  for (vector<StochasticNode*>::const_iterator p(_nodes.begin());
       p != _nodes.end(); ++p) 
    {
      double l = (*p)->logDensity(_chain);
      if (isnan(l)) {
	throw NodeError(*p, "Failure to calculate log density");
      }
      else if (l == -DBL_MAX || l == DBL_MAX || isinf(l)) {
        return l;
      }
      else {
	logprior += l;
      }
    }

  double loglikelihood = 0;
  for (vector<StochasticNode const*>::const_iterator p(_stoch_children.begin());
       p != _stoch_children.end(); ++p) 
    {
      double l = (*p)->logDensity(_chain);
      if (isnan(l)) {
	  throw NodeError(*p, "Failure to calculate log density");
      }
      else if (l == -DBL_MAX || l == DBL_MAX || isinf(l)) {
	return l;
      }
      else {
	loglikelihood += l;
      }
    }

  return logprior + loglikelihood;
}

vector<StochasticNode const*> const &Sampler::stochasticChildren() const
{
  return _stoch_children;
}

vector<DeterministicNode*> const &Sampler::deterministicChildren() const
{
  return _determ_children;
}

void Sampler::setValue(vector<double const *> const &value, 
		       vector<unsigned long> const &length)
{
    unsigned int n = _nodes.size();
    if (value.size() != n || length.size() != n) {
      throw logic_error("Argument length mismatch in Sampler::setValue");
    }
    for (unsigned int i = 0; i < n; ++i) {
      _nodes[i]->setValue(value[i], length[i], _chain);
    }
    for (vector<DeterministicNode*>::iterator p(_determ_children.begin());
	 p != _determ_children.end(); ++p) {
      (*p)->deterministicSample(_chain);
    }
}

unsigned int Sampler::chain()
{
  return _chain;
}
