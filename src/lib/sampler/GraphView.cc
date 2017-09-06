#include <config.h>
#include <sampler/GraphView.h>
#include <graph/StochasticNode.h>
#include <graph/DeterministicNode.h>
#include <graph/Graph.h>
#include <graph/NodeError.h>
#include <util/nainf.h>

#include <stdexcept>
#include <set>
#include <list>
#include <string>
#include <cmath>
#include <algorithm>

using std::vector;
using std::set;
using std::list;
using std::runtime_error;
using std::logic_error;
using std::string;
using std::copy;

static unsigned int sumLength(vector<jags::StochasticNode *> const &nodes)
{
    //Adds up the length of a vector of stochastic nodes
    unsigned int n = 0;
    for (unsigned int i = 0; i < nodes.size(); ++i) {
	n += nodes[i]->length();
    }
    return n;
}

namespace jags {

GraphView::GraphView(vector<StochasticNode *> const &nodes, Graph const &graph,
		     bool multilevel)
    : _length(sumLength(nodes)), _nodes(nodes), _stoch_children(0),
      _determ_children(0), _multilevel(false)
{
    //Sanity check on node
    //FIXME: Could use a templated version of countChains here
    for (unsigned int i = 1; i < nodes.size(); ++i) {
	if (nodes[i]->nchain() != nodes[0]->nchain()) {
	    throw logic_error("Chain mismatch in GraphView");
	}
    }
    classifyChildren(nodes, graph, _stoch_children, _determ_children,
		     multilevel);
}

vector<StochasticNode *> const &GraphView::nodes() const
{
  return _nodes;
}

static bool classifyNode(StochasticNode *snode, Graph const &sample_graph, 
			 set<StochasticNode const *> &sset,
			 list<StochasticNode *> &slist)
{
    // classification function for stochastic nodes

    if (sset.count(snode))
	return true;
    
    if (sample_graph.contains(snode)) {
	sset.insert(snode);
	slist.push_back(snode);
	return true;
    }
    else {
	return false;
    }
}

static bool classifyNode(DeterministicNode *dnode, 
			 Graph const &sample_graph,
                         set<StochasticNode const *> &sset,
			 list<StochasticNode *> &slist,
			 set<DeterministicNode const *> &dset,
			 list<DeterministicNode *> &dlist)
{
    //  Recursive classification function for deterministic nodes

    if (!sample_graph.contains(dnode))
	return false;
    
    if (dset.count(dnode))
	return true;
    
    bool informative = false;
    list<StochasticNode*>::const_iterator p; 
    for (p = dnode->stochasticChildren()->begin(); 
	 p != dnode->stochasticChildren()->end(); ++p)
    {
	if (classifyNode(*p, sample_graph, sset, slist))
	    informative = true;
    }
    list<DeterministicNode*>::const_iterator q;
    for (q = dnode->deterministicChildren()->begin();
	 q != dnode->deterministicChildren()->end(); ++q)
    {
	if (classifyNode(*q, sample_graph, sset, slist, dset, dlist)) 
	    informative = true;
    }
    if (informative) {
	dset.insert(dnode);
	dlist.push_back(dnode);
    }
    return informative;
}


void GraphView::classifyChildren(vector<StochasticNode *> const &nodes,
				 Graph const &graph,
				 vector<StochasticNode *> &stoch_nodes,
				 vector<DeterministicNode*> &dtrm_nodes,
				 bool multilevel)
{
    set<StochasticNode const *> sset;
    set<DeterministicNode const *> dset;
    list<StochasticNode *> slist;
    list<DeterministicNode *> dlist;

    /* Classify children of each node */
    vector<StochasticNode  *>::const_iterator p; 
    for (p = nodes.begin(); p != nodes.end(); ++p) {
	if (!graph.contains(*p)) {
	    throw logic_error("Sampled node outside of sampling graph");
	}
	list<StochasticNode*> const *sch = (*p)->stochasticChildren();
	for (list<StochasticNode*>::const_iterator q = sch->begin();
	     q != sch->end(); ++q)
	{
	    classifyNode(*q, graph, sset, slist);
	}
	list<DeterministicNode*> const *dch = (*p)->deterministicChildren();
	for (list<DeterministicNode*>::const_iterator q = dch->begin();
	     q != dch->end(); ++q)
	{
	    classifyNode(*q, graph, sset, slist, dset, dlist);
	}
    }

    if (multilevel) {
	/* Strip nodes to be sampled out of the set of stochastic
	   children. Such nodes would contribute to both the prior
	   AND the likelihood, causing incorrect calculation of the
	   log full conditional */
	for (p = nodes.begin(); p != nodes.end(); ++p) {
	    if (sset.count(*p)) {
		list<StochasticNode*>::iterator i = 
		    find(slist.begin(), slist.end(), *p);
		if (i == slist.end()) {
		    throw logic_error("error in ClassifyChildren"); 
		}
		else {
		    slist.erase(i);
		}
	    }
	}
	/* 
	   We also need ensure that we calculate the full log density
	   for each sampled node.
	*/
	_multilevel = true;
    }
    else {
	for (p = nodes.begin(); p != nodes.end(); ++p) {
	    if (sset.count(*p)) {
		throw logic_error("Invalid multilevel GraphView");
	    }
	}
	
    }

    stoch_nodes.clear();
    for (list<StochasticNode *>::const_iterator i = slist.begin();
         i != slist.end(); ++i)
    {
	stoch_nodes.push_back(*i);
    }

    dtrm_nodes.clear();
    // Deterministic nodes are pushed onto dtrm_nodes in reverse order
    for (list<DeterministicNode *>::reverse_iterator i = dlist.rbegin();
         i != dlist.rend(); ++i)
    {
	dtrm_nodes.push_back(*i);
    }

}

double GraphView::logFullConditional(unsigned int chain) const
{
    PDFType pdf_prior = _multilevel ? PDF_FULL : PDF_PRIOR;

    double lprior = 0.0;
    vector<StochasticNode*>::const_iterator p = _nodes.begin();
    for (; p != _nodes.end(); ++p) {
	lprior += (*p)->logDensity(chain, pdf_prior);
    }
  
    double llike = 0.0;
    vector<StochasticNode *>::const_iterator q = _stoch_children.begin();
    for (; q != _stoch_children.end(); ++q) {
	llike += (*q)->logDensity(chain, PDF_LIKELIHOOD);
    }

    double lfc = lprior + llike;
    if(jags_isnan(lfc)) {
	/* 
	   Try to find where the calculation went wrong. At this point,
	   we are definitely going to throw an exception. It's just a
	   question of working out which error message. So we can afford
	   to be laborious.
	*/

	//Check prior
	for (p = _nodes.begin(); p != _nodes.end(); ++p) {
	    if (jags_isnan((*p)->logDensity(chain, pdf_prior))) {
		throw NodeError(*p, "Failure to calculate log density");
	    }
	}
	if (jags_isnan(lprior)) {
	    throw runtime_error("Failure to calculate prior density");
	}

	//Recalculate the deterministic children, checking for
	//invalid values
	vector<DeterministicNode*>::const_iterator r;
	for (r =_determ_children.begin(); r != _determ_children.end(); ++r) {
	    if(!(*r)->checkParentValues(chain)) {
		throw NodeError(*r, "Invalid parent values");
	    }
	    (*r)->deterministicSample(chain);
	}

	//Check likelihood
	for (q = _stoch_children.begin(); q != _stoch_children.end(); ++q) {
	    if (jags_isnan((*q)->logDensity(chain, PDF_LIKELIHOOD))) {
		throw NodeError(*q, "Failure to calculate log density");
	    }
	}
	if (jags_isnan(llike)) {
	    throw runtime_error("Failure to calculate likelihood");
	}

	//This could happen adding -Inf to +Inf
	if (!jags_finite(lprior) && !jags_finite(llike)) {
	    throw runtime_error("Prior and likelihood are incompatible");
	}

	//Something else went wrong, but what?
	throw runtime_error("Failure to calculate log full conditional");
    }

    return lfc;
}

double GraphView::logPrior(unsigned int chain) const
{
    //In a multi-level GraphView we need to calculate the full log
    //density of each sampled node
    PDFType pdf_prior = _multilevel ? PDF_FULL : PDF_PRIOR;

    double lprior = 0.0;
    vector<StochasticNode*>::const_iterator p = _nodes.begin();
    for (; p != _nodes.end(); ++p) {
	lprior += (*p)->logDensity(chain, pdf_prior);
    }
  
    if(jags_isnan(lprior)) {
	//Try to find where the calculation went wrong
	for (p = _nodes.begin(); p != _nodes.end(); ++p) {
	    if (jags_isnan((*p)->logDensity(chain, pdf_prior))) {
		throw NodeError(*p, "Failure to calculate log prior density");
	    }
	}
	throw logic_error("Failure in GraphView::logLikelihood");
    }

    return lprior;
}

double GraphView::logLikelihood(unsigned int chain) const
{
    double llik = 0.0;

    vector<StochasticNode *>::const_iterator q = _stoch_children.begin();
    for (; q != _stoch_children.end(); ++q) {
	llik += (*q)->logDensity(chain, PDF_LIKELIHOOD);
    }
  
    if(jags_isnan(llik)) {
	//Try to find where the calculation went wrong
	for (q = _stoch_children.begin(); q != _stoch_children.end(); ++q) {
	    if (jags_isnan((*q)->logDensity(chain, PDF_LIKELIHOOD))) {
		throw NodeError(*q, "Failure to calculate log likelihood");
	    }
	}

	//This could  happen if we try to add +Inf to -Inf
	throw logic_error("Failure in GraphView::logLikelihood");
    }

    return llik;
}

vector<StochasticNode *> const &GraphView::stochasticChildren() const
{
  return _stoch_children;
}

vector<DeterministicNode*> const &GraphView::deterministicChildren() const
{
  return _determ_children;
}

void GraphView::setValue(double const * value, unsigned int length,
			 unsigned int chain) const
{
    if (length != _length) {
      throw logic_error("Argument length mismatch in GraphView::setValue");
    }

    for (unsigned int i = 0; i < _nodes.size(); ++i) {
	Node *node = _nodes[i];
	node->setValue(value, node->length(), chain);
	value += node->length();
    }

    for (vector<DeterministicNode*>::const_iterator p(_determ_children.begin());
	 p != _determ_children.end(); ++p) {
      (*p)->deterministicSample(chain);
    }
}

void GraphView::setValue(vector<double> const &value, unsigned int chain) const
{
    //FIXME We could inline this
    setValue(&value[0], value.size(), chain);
}
 
void GraphView::getValue(vector<double> &value, unsigned int chain) const 
{
    if (value.size() != _length) 
	throw logic_error("length mismatch in GraphView::getValue");
    
    unsigned int k = 0;
    for (unsigned int i = 0; i < _nodes.size(); ++i) {
	double const *vi = _nodes[i]->value(chain);
	for (unsigned int j = 0; j < _nodes[i]->length(); ++j) {
	    value[k++] = vi[j];
	}
    }
}

unsigned int GraphView::length() const
{
    return _length;

}

bool GraphView::isDependent(Node const *node) const
{
    for (unsigned int i = 0; i < _nodes.size(); ++i) {
	if (_nodes[i] == node)
	    return true;
    }
    for (unsigned int j = 0; j < _determ_children.size(); ++j) {
	if (_determ_children[j] == node)
	    return true;
    }
    return false;
}
      
unsigned int nchain(GraphView const *gv)
{
    return gv->nodes()[0]->nchain();
}

    
    void GraphView::checkFinite(unsigned int chain) const
    {

	vector<StochasticNode*>::const_iterator p = _nodes.begin();
	for ( ; p != _nodes.end(); ++p)
	{
	    double ld = (*p)->logDensity(chain, PDF_PRIOR);
	    if (jags_isnan(ld)) {
		throw NodeError(*p, "Error calculating log density");
	    }
	    else if (ld == JAGS_NEGINF || (!jags_finite(ld) && ld < 0)) {
		throw NodeError(*p, "The value of this node is inconsistent with its prior distribution.");
	    }
	}

	
	for (p =_stoch_children.begin(); p != _stoch_children.end(); ++p)
	{
	    double ld = (*p)->logDensity(chain, PDF_PRIOR);
	    if (jags_isnan(ld)) {
		throw NodeError(*p, "Error calculting log density");
	    }
	    else if (ld == JAGS_NEGINF || (!jags_finite(ld) && ld < 0)) {
		throw NodeError(*p, "Node inconsistent with parents");
	    }
	}

    }

} //namespace jags
