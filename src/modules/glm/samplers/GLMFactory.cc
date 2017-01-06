#include <config.h>

#include "GLMFactory.h"
#include "GLMSampler.h"

#include <graph/GraphMarks.h>
#include <graph/Graph.h>
#include <graph/StochasticNode.h>
#include <graph/DeterministicNode.h>
#include <graph/LinkNode.h>
#include <distribution/Distribution.h>
#include <sampler/Linear.h>
#include <sampler/GraphView.h>
#include <sampler/SingletonGraphView.h>

#include <set>
#include <map>
#include <algorithm>
#include <utility>

using std::set;
using std::vector;
using std::map;
using std::stable_sort;
using std::string;
using std::list;
using std::pair;

namespace jags {

/*
  Aggregates candidate Nodes with common stochastic children. 
*/
static bool aggregate(SingletonGraphView const *candidate,
		      set<StochasticNode const *> &stochastic_children,
		      Graph const &graph)
{

    // Check that there is some overlap in stochastic children between
    // candidate node and current set.  

    vector<StochasticNode *> const &candidate_children = 
	candidate->stochasticChildren();

    bool overlap = false;
    for (unsigned int i = 0; i < candidate_children.size(); ++i) {
	if (stochastic_children.count(candidate_children[i]) > 0) {
	    stochastic_children.insert(candidate_children.begin(), 
				       candidate_children.end());
	    return true;
	}
    }
    
    return overlap;
}

    typedef pair<SingletonGraphView*,unsigned int> ViewScore;
    
    struct less_viewscore {
	/* 
	   Comparison operator for ViewScore objects that ranks them
	   in order of the score (second member of the pair);
	*/
	
	bool operator()(ViewScore const &x, ViewScore const &y) const {
	    return x.second < y.second;
	};
    };
    
    vector<SingletonGraphView*>
    sort_by_children(vector<SingletonGraphView*> const &candidates)
    {
	// Utility function called by makeSampler when we have a set
	// of candidates for sampling that are not jointly linear. To
	// find a maximal subset of linear candidates we first sort in
	// increasing order of the number of shared outcome
	// variables. The heuristic reasoning is that a candidate with
	// few shared outcomes has few opportunities for a quadratic
	// interaction with another candidate.

	map<StochasticNode*,unsigned int> cmap; // Score map for candidates

	// First pass. Calculate scores for stochastic children equal
	// to the number of candidates that they depend on, minus 1. 
	for (vector<SingletonGraphView*>::const_iterator p = candidates.begin();
	     p != candidates.end(); ++p)
	{
	    vector<StochasticNode*> const &sc = (*p)->stochasticChildren();
	    for (vector<StochasticNode*>::const_iterator q = sc.begin();
		 q != sc.end(); ++q)
	    {
		map<StochasticNode*,unsigned int>::iterator s = cmap.find(*q);
		if (s == cmap.end()) {
		    cmap[*q] = 0;
		}
		else {
		    s->second++;
		}
	    }
	}

	vector<ViewScore> vscores;

	// Second pass. For each candidate, add up the scores of
	// stochastic children. This is the number of stochastic
	// children shared with other candidates.
	for (vector<SingletonGraphView*>::const_iterator p = candidates.begin();
	     p != candidates.end(); ++p)
	{
	    ViewScore pshared(*p,0);
	    vector<StochasticNode*> const &sc = (*p)->stochasticChildren();
	    for (vector<StochasticNode*>::const_iterator q = sc.begin();
		 q != sc.end(); ++q)
	    {
		map<StochasticNode*,unsigned int>::iterator s = cmap.find(*q);
		pshared.second += s->second;
	    }
	    vscores.push_back(pshared);
	}

	// Now sort candidates in score order
	stable_sort(vscores.begin(), vscores.end(), less_viewscore());
	vector<SingletonGraphView*> sorted_candidates;
	for (vector<ViewScore>::const_iterator p = vscores.begin();
	     p != vscores.end(); ++p)
	{
	    sorted_candidates.push_back(p->first);
	}
	return sorted_candidates;
    }

    static bool isLink(DeterministicNode const *dnode)
    {
	return dynamic_cast<LinkNode const*>(dnode) != 0;
    }
    
    vector<bool>
    get_linear_subset(vector<SingletonGraphView*> const &candidates,
		      bool fixed, bool link)
    {
	vector<bool> linear;
	
	set<Node const *> ancestors;
	set<Node const *> link_children;
	pair<set<Node const*>::iterator, bool> p;
	
	for (vector<SingletonGraphView*>::const_iterator c = candidates.begin();
	     c != candidates.end(); ++c)
	{
	    vector <set<Node const*>::iterator> new_ancestors;
	    vector <set<Node const*>::iterator> new_link_children;
	    
	    StochasticNode *snode = (*c)->node();
	    p = ancestors.insert(snode);
	    new_ancestors.push_back(p.first);

	    bool ok = true;
	    vector<DeterministicNode*> const &dn =
		(*c)->deterministicChildren();

	    for (unsigned int j = 0; j < dn.size(); ++j) {

		if (link_children.count(dn[j])) {
		    ok = false;
		    break;
		}
		
		if (dn[j]->isClosed(ancestors, DNODE_LINEAR, fixed)) {
		    p = ancestors.insert(dn[j]);
		    if (p.second) {
			new_ancestors.push_back(p.first);
		    }
		}
		else if (link && isLink(dn[j])) {
		    list<DeterministicNode*>::const_iterator i;
		    for (i = dn[j]->deterministicChildren()->begin();
			 i != dn[j]->deterministicChildren()->end(); ++i)
		    {
			p = link_children.insert(*i);
			if (p.second) {
			    new_link_children.push_back(p.first);
			}
		    }
		}
		else {
		    ok = false;
		    break;
		}
	    }

	    if (!ok) {

		// Non-linear interaction found. Unwind additions to
		// ancestors and link_children
		
		for (unsigned int i = 0; i < new_ancestors.size(); ++i) {
		    ancestors.erase(new_ancestors[i]);
		}
		for (unsigned int i = 0; i < new_link_children.size(); ++i) {
		    link_children.erase(new_link_children[i]);
		}
	    }
	    linear.push_back(ok);
	}

	return linear;
    }


namespace glm {
    
    GLMFactory::GLMFactory(string const &name)
	: _name(name)
    {}

    bool GLMFactory::checkDescendants(SingletonGraphView const *view) const
    {
	// Check stochastic children
	vector<StochasticNode *> const &stoch_nodes = 
	    view->stochasticChildren();
	for (unsigned int i = 0; i < stoch_nodes.size(); ++i) {
	    if (isBounded(stoch_nodes[i])) {
		return false; //Truncated outcome variable
	    }
	    if (!checkOutcome(stoch_nodes[i])) {
		return false; //Invalid outcome or link
	    }
	    if (fixedOutcome() && !isObserved(stoch_nodes[i])) {
		return false; //Unobserved outcome not allowed by sampler
	    }
	    //Check that other parameters do not depend on snode	    
	    vector<Node const *> const &param = stoch_nodes[i]->parents();
	    for (unsigned int j = 1; j < param.size(); ++j) {
		if (view->isDependent(param[j])) {
		    return false;
		}
	    }
	}

	// Check linearity of deterministic descendants
	if (!checkLinear(view, fixedDesign(), true))
	    return false;

	return true;
    }


    SingletonGraphView * 
    GLMFactory::makeView(StochasticNode *snode, Graph const &graph,
			 bool gibbs) const
    {
	/*
	  Returns a newly allocated GraphView if node can be sampled,
	  otherwise zero pointer.
	*/

	string dname = snode->distribution()->name();
	if (dname != "dnorm" && dname != "dmnorm")
	    return 0; //Must have normal prior

	if (gibbs && snode->length() != 1) {
	    return 0;
	}
	if (!gibbs && isBounded(snode)) {
	    return 0;
	}

	SingletonGraphView *view = new SingletonGraphView(snode, graph);
	if (!checkDescendants(view)) {
	    delete view;
	    return 0;
	}
	else {
	    return view;
	}
    }

    GLMFactory::~GLMFactory()
    {}
    
    Sampler * 
    GLMFactory::makeSampler(list<StochasticNode*> const &free_nodes, 
			    Graph const &graph, bool gibbs) const
    {
	// Find candidate nodes that could be part of a GLM.
	vector<SingletonGraphView*> candidates;
	for (list<StochasticNode*>::const_iterator p = free_nodes.begin();
	     p != free_nodes.end(); ++p)
	{
	    if (SingletonGraphView *up = makeView(*p, graph, gibbs)) {
		candidates.push_back(up);
	    }
	}
	if (candidates.empty()) {
	    return 0;
	}

	// Now try to aggregate nodes into a joint linear model
	unsigned int Nc = candidates.size();
	vector<bool> keep(Nc, false);
	vector<bool> resolved(Nc, false);
	GraphView *view = 0;
	    
	keep[0] = true;
	
	set<StochasticNode const *> stoch_children;
	stoch_children.insert(candidates[0]->stochasticChildren().begin(),
			      candidates[0]->stochasticChildren().end());

	// Aggregate nodes with common stochastic children.
	// This requires multiple passes through candidate nodes
	bool loop = true;
	while (loop) {
	    loop = false;
	    for (unsigned int j = 1; j < candidates.size(); ++j) {
		if (!keep[j]) {
		    keep[j] = aggregate(candidates[j], stoch_children, graph);
		    if (keep[j]) {
			loop = true;
		    }
		}
	    }
	}

	// Remove candidate nodes that are stochastic children of
	// another candidate node. All GLMMethod algorithms rely on
	// the stochastic parents of the sampled nodes being fixed
	// within any update.

	set<StochasticNode const*> all_children;
	for (unsigned int j = 0; j < candidates.size(); ++j) {
	    if (keep[j]) {
		vector<StochasticNode *> const &children_j =
		    candidates[j]->stochasticChildren();
		all_children.insert(children_j.begin(), children_j.end());
	    }
	}
	for (unsigned int j = 0; j < candidates.size(); ++j) {
	    if (keep[j]) {
		if (all_children.count(candidates[j]->node())) {
		    keep[j] = false;
		}
	    }
	}

	//
	vector<StochasticNode*> sample_nodes;
	vector<SingletonGraphView*> sub_views;
	for (unsigned int j = 0; j < Nc; ++j) {
	    if (keep[j]) {
		sub_views.push_back(candidates[j]);
		sample_nodes.push_back(candidates[j]->node());
	    }
	    else {
		delete candidates[j];
	    }
	}

	if (sub_views.empty()) {
	    return 0;
	}

	view = new GraphView(sample_nodes, graph);
	if (!checkLinear(view, fixedDesign(), true)) {

	    // Candidates are not jointly linear
	    
	    sample_nodes.clear();
	    sub_views.clear();
	    delete view; 

	    // Find maximal linear subset
	    
	    vector<SingletonGraphView*> sorted_candidates = 
		sort_by_children(candidates);
	    vector<bool> keep = get_linear_subset(sorted_candidates,
						  fixedDesign(), true);
	    for (unsigned int j = 0; j < sorted_candidates.size(); ++j) {
		if (keep[j]) {
		    sub_views.push_back(sorted_candidates[j]);
		    sample_nodes.push_back(sorted_candidates[j]->node());
		}
		else {
		    delete sorted_candidates[j];
		}
	    }
	    view = new GraphView(sample_nodes, graph);
	}
	
	/* FIXME
	 * 
	 * The overhead of a GLM sampling method for a univariate
	 * sampler is too large So we skip univariate nodes. However,
	 * by returning a null pointer at this point we terminate the
	 * factory and so will miss a block of nodes that could be
	 * sampled together later.
	 *
	 * What we need to do here is to mark this node to be skipped in any
	 * later call.
	 */
	if (sub_views.size() == 1 && view->length() == 1) {
	    delete view;
	    delete sub_views.back();
	    return 0;
	}
	
	unsigned int Nch = nchain(view);
	vector<MutableSampleMethod*> methods(Nch, 0);
	
	vector<SingletonGraphView const*> const_sub_views(sub_views.size());
	copy(sub_views.begin(), sub_views.end(), const_sub_views.begin());
	for (unsigned int ch = 0; ch < Nch; ++ch) {
	    methods[ch] = newMethod(view, const_sub_views, ch, gibbs);
	}

	string nm = _name;
	if (gibbs) nm.append("-Gibbs");
	return new GLMSampler(view, sub_views, methods, nm);
    }

    string GLMFactory::name() const
    {
	return _name;
    }

    vector<Sampler*>  
    GLMFactory::makeSamplers(list<StochasticNode*> const &nodes, 
			     Graph const &graph) const
    {
	if (Sampler *s = makeSampler(nodes, graph, false)) {
	    return vector<Sampler*>(1, s);
	}
	if (Sampler *s = makeSampler(nodes, graph, true)) {
	    return vector<Sampler*>(1, s);
	}
	return vector<Sampler*>();
    }

    bool GLMFactory::fixedDesign() const
    {
	return false;
    }

    bool GLMFactory::fixedOutcome() const
    {
	return false;
    }
    
}}

