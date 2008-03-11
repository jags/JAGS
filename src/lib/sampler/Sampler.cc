#include <config.h>
#include <sampler/Sampler.h>
#include <graph/StochasticNode.h>
#include <graph/Graph.h>
#include <graph/NodeError.h>
#include <util/nainf.h>

#include <stdexcept>
#include <set>
#include <string>
#include <cmath>

using std::vector;
using std::set;
using std::runtime_error;
using std::logic_error;
using std::string;

static unsigned int sumLength(vector<StochasticNode *> const &nodes)
{
    //Adds up the length of a vector of stochastic nodes
    unsigned int n = 0;
    for (unsigned int i = 0; i < nodes.size(); ++i) {
	n += nodes[i]->length();
    }
    return n;
}

Sampler::Sampler(vector<StochasticNode *> const &nodes, Graph const &graph)
    : _length(sumLength(nodes)), _nodes(nodes), _stoch_children(0),
      _determ_children(0)
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

static bool isInformative(Node *node, Graph const &sample_graph)
{
    /* Recursively test whether a node is informative (is observed or
       has an observed descendant in the graph) */

    if (!sample_graph.contains(node))
	return false;
    
    if (node->isObserved())
	return true;

    for (set<Node*>::iterator p = node->children()->begin();
	 p != node->children()->end(); ++p) {
	if (isInformative(*p, sample_graph)) {
	    return true;
	}
    }    
    return false;
}

static bool classifyNode(Node *node, Graph const &sample_graph, 
			 Graph &sgraph, Graph &dgraph)
{
    /*
     * Recursive classification function for node and its descendants.
     *
     * If a node is informative, it is added to either sgraph (for
     * stochastic nodes representing random variables) or dgraph (
     * for other nodes). The recursion stops with stochastic nodes.
     */
    if (!sample_graph.contains(node))
	return false;

    bool isinformative = false;
    if (asStochastic(node) && node->isRandomVariable()) {
	/* This might look redundant, but we want to leave open
	   the possibility that stochastic nodes may not
	   represent random variables */

	if (sgraph.contains(node)) 
	    return true;

	isinformative = isInformative(node, sample_graph);
	if (isinformative) {
	    sgraph.add(node);
	}
    }
    else {

	if (dgraph.contains(node))
	    return true;
	
	for (set<Node*>::iterator p = node->children()->begin();
	     p != node->children()->end(); ++p) {
	    if (classifyNode(*p, sample_graph, sgraph, dgraph)) {
		isinformative = true;
	    }
	}
	if (isinformative) {
	    dgraph.add(node);
	}
    }
    return isinformative;
}

void Sampler::classifyChildren(vector<StochasticNode *> const &nodes,
			       Graph const &graph,
			       vector<StochasticNode const*> &stoch_nodes,
			       vector<Node*> &dtrm_nodes)
{
    Graph dgraph, sgraph;

    /* Classify children of each node */
    vector<StochasticNode  *>::const_iterator p = nodes.begin();
    for (; p != nodes.end(); ++p) {
	StochasticNode *snode = *p;
	if (!graph.contains(snode)) {
	    throw logic_error("Sampled node outside of sampling graph");
	}
	for (set<Node*>::const_iterator q = snode->children()->begin(); 
	     q != snode->children()->end(); ++q) 
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
    stoch_nodes.clear();
    for (vector<Node*>::iterator i = svector.begin(); i != svector.end(); 
	 ++i)
    {
	stoch_nodes.push_back(asStochastic(*i));
    }

    dtrm_nodes.clear();
    dgraph.getSortedNodes(dtrm_nodes);
}

double Sampler::logFullConditional(unsigned int chain) const
{
    double lprior = 0.0;
    vector<StochasticNode*>::const_iterator p = _nodes.begin();
    for (; p != _nodes.end(); ++p) {
	lprior += (*p)->logDensity(chain);
    }
  
    double llike = 0.0;
    vector<StochasticNode const*>::const_iterator q = _stoch_children.begin();
    for (; q != _stoch_children.end(); ++q) {
	llike += (*q)->logDensity(chain);
    }

    double lfc = lprior + llike;
    if(jags_isnan(lfc)) {
	/* 
	   Try to find where the calculation went wrong. At this point,
	   we are definitely going to throw an exception. It's just a
	   question of working out which error message. So we can afford
	   to be laborious.
	*/
	for (p = _nodes.begin(); p != _nodes.end(); ++p) {
	    if (jags_isnan((*p)->logDensity(chain))) {
		throw NodeError(*p, "Failure to calculate log density");
	    }
	}
    
	for (q = _stoch_children.begin(); q != _stoch_children.end(); ++q) {
	    if (jags_isnan((*q)->logDensity(chain))) {
		throw NodeError(*q, "Failure to calculate log density");
	    }
	}

	if (jags_isnan(lprior)) {
	    throw runtime_error(string("Failure to calculate prior density in ")
				+ name());
	}
	if (jags_isnan(llike)) {
	    throw runtime_error(string("Failure to calculate likelihood in ")
				+ name());
	}
	if (!jags_finite(lprior) && !jags_finite(llike)) {
	    throw runtime_error(string("Prior and likelihood are incompatible")
				+ " in " + name());
	}
	throw runtime_error(string("Failure to calculate log full conditional")
			    + " in " + name());
    }

    return lfc;
}

double Sampler::logPrior(unsigned int chain) const
{
    double lprior = 0.0;

    vector<StochasticNode*>::const_iterator p = _nodes.begin();
    for (; p != _nodes.end(); ++p) {
	lprior += (*p)->logDensity(chain);
    }
  
    if(jags_isnan(lprior)) {
	//Try to find where the calculation went wrong
	for (p = _nodes.begin(); p != _nodes.end(); ++p) {
	    if (jags_isnan((*p)->logDensity(chain))) {
		throw NodeError(*p, "Failure to calculate log density");
	    }
	}
	throw logic_error("Failure in Sampler::logLikelihood");
    }

    return lprior;
}

double Sampler::logLikelihood(unsigned int chain) const
{
    double llik = 0.0;

    vector<StochasticNode const*>::const_iterator q = _stoch_children.begin();
    for (; q != _stoch_children.end(); ++q) {
	llik += (*q)->logDensity(chain);
    }
  
    if(jags_isnan(llik)) {
	//Try to find where the calculation went wrong
	for (q = _stoch_children.begin(); q != _stoch_children.end(); ++q) {
	    if (jags_isnan((*q)->logDensity(chain))) {
		throw NodeError(*q, "Failure to calculate log density");
	    }
	}

	//This could  happen if we try to add +Inf to -Inf
	throw logic_error("Failure in Sampler::logLikelihood");
    }

    return llik;
}

vector<StochasticNode const*> const &Sampler::stochasticChildren() const
{
  return _stoch_children;
}

vector<Node*> const &Sampler::deterministicChildren() const
{
  return _determ_children;
}

void Sampler::setValue(double const * value, unsigned int length,
		       unsigned int chain)
{
    if (length != _length) {
      throw logic_error("Argument length mismatch in Sampler::setValue");
    }

    for (unsigned int i = 0; i < _nodes.size(); ++i) {
	Node *node = _nodes[i];
	node->setValue(value, node->length(), chain);
	value += node->length();
    }

    for (vector<Node*>::iterator p(_determ_children.begin());
	 p != _determ_children.end(); ++p) {
      (*p)->deterministicSample(chain);
    }
}

unsigned int Sampler::length() const
{
    return _length;
}
