#include <config.h>
#include <sampler/Linear.h>
#include <sampler/GraphView.h>
#include <sampler/SingletonGraphView.h>
#include <graph/GraphMarks.h>
#include <graph/Graph.h>
#include <graph/StochasticNode.h>
#include <graph/DeterministicNode.h>
#include <graph/LinkNode.h>
#include <graph/Node.h>

using std::vector;
using std::set;
using std::list;

namespace jags {

typedef bool (DeterministicNode::*NodeCheckFn) (GraphMarks const&, bool) const;

static bool isLink(DeterministicNode const *dnode)
{
    return dynamic_cast<LinkNode const*>(dnode) != nullptr;
}

    bool checkAdditive(SingletonGraphView const *sgv, bool fixed)
    {
	vector<DeterministicNode *> const &dn =	sgv->deterministicChildren();
	set<Node const*> ancestors;
	ancestors.insert(sgv->node());
	for (unsigned int i = 0; i < dn.size(); ++i) {
	    if(dn[i]->isClosed(ancestors, DNODE_ADDITIVE, fixed)) {
		ancestors.insert(dn[i]);
	    }
	    else {
		return false;
	    }
	}
	return true;
    }

    bool checkAdditive(vector<StochasticNode*> const &snodes,
		       Graph const &graph, bool fixed)
    {
	//A GraphView is additive if every SingletonGraphview for each
	//sampled node in the same graph is additive.
	for (unsigned int i = 0; i < snodes.size(); ++i) {
	    SingletonGraphView sgv(snodes[i], graph);
	    if (!checkAdditive(&sgv, false)) return false;
	}

	if (!fixed) return true;

	//If we want to ensure that the additive function has a fixed
	//intercept, we need to ensure that all parent nodes of
	//deterministic children are either additive nodes themselves
	//or are fixed.
	
	set<Node const*> ancestors;
	//Sampled nodes are trivial (fixed) additive functions of themselves
#ifndef _RWSTD_NO_MEMBER_TEMPLATES
	ancestors.insert(snodes.begin(), snodes.end());
#else
	//Workaround for Solaris libCstd
	for (vector<StochasticNode*>::const_iterator p =
		 snodes.begin(); p != snodes.end(); ++p)
	{
	    ancestors.insert(*p);
	}
#endif

	GraphView gv(snodes, graph);
	vector<DeterministicNode *> const &dn =	gv.deterministicChildren();
	for (unsigned int j = 0; j < dn.size(); ++j) {
	    vector<Node const *> parj = dn[j]->parents();
	    for (unsigned int k = 0; k < parj.size(); ++k) {
		if (ancestors.count(parj[k]) == 0 && !parj[k]->isFixed()) {
		    return false;
		}
	    }
	    ancestors.insert(dn[j]);
	}
	return true;
    }
    
bool checkLinear(GraphView const *gv, bool fixed, bool link)
{
    vector<DeterministicNode *> const &dn = gv->deterministicChildren();

    set<Node const*> ancestors;
    //Sampled nodes are trivial (fixed) linear functions of themselves
#ifndef _RWSTD_NO_MEMBER_TEMPLATES
    ancestors.insert(gv->nodes().begin(), gv->nodes().end());
#else
    //Workaround for Solaris libCstd
    for (vector<StochasticNode const *>::const_iterator p = gv->nodes().begin();
         p != gv->nodes().end(); ++p)
      {
        ancestors.insert(*p);
      }
#endif
    
    for (unsigned int j = 0; j < dn.size(); ++j) {
	if (dn[j]->isClosed(ancestors, DNODE_LINEAR, fixed)) {
	    ancestors.insert(dn[j]);
	}
	else if (link && isLink(dn[j])) {
	    // A link function is allowed if no other deterministic
	    // nodes in the GraphView depend on it.
	    set<DeterministicNode*> dset;
	    list<DeterministicNode*> const *dc = dn[j]->deterministicChildren();
#ifndef _RWSTD_NO_MEMBER_TEMPLATES
	    dset.insert(dc->begin(), dc->end());
#else
	    for (list<DeterministicNode*>::const_iterator p = dc->begin(); 
		 p != dc->end(); ++p) 
            {
	        dset.insert(*p);
	    }
#endif
	    for (unsigned int k = j + 1; k < dn.size(); ++k) {
		if (dset.count(dn[k])) {
		    return false;
		}
	    }
	}
	else {
	    return false;
	}
    }
    
    return true;
}

bool checkScale(GraphView const *gv, bool fixed)
{
    vector<DeterministicNode *> const &dnodes = gv->deterministicChildren();
    
    set<Node const*> ancestors;
    //FIXME: valid for multiple nodes?
#ifdef _RWSTD_NO_MEMBER_TEMPLATES
    //Workaround for Solaris libCstd
    for (vector<StochasticNode const *>::const_iterator p = gv->nodes().begin();
         p != gv->nodes().end(); ++p)
      {
        ancestors.insert(*p);
      }
#else
    ancestors.insert(gv->nodes().begin(), gv->nodes().end());
#endif 

    //Start off looking for scale transformations, then fall back on
    //scale mixture transformations if fixed is false.

    bool mix = false;
    for (unsigned int j = 0; j < dnodes.size(); ++j) {

	if (!mix) {
	    if (dnodes[j]->isClosed(ancestors, DNODE_SCALE, fixed)) {
		ancestors.insert(dnodes[j]);
	    }
	    else if (!fixed) {
		mix = true;
	    }
	    else {
		return false;
	    }
	}

	if (mix) {
	    if (dnodes[j]->isClosed(ancestors, DNODE_SCALE_MIX, false)) {
		ancestors.insert(dnodes[j]);
	    }
	    else {
		return false;
	    }
	}
    }
    
    return true;
}

bool checkPower(GraphView const *gv, bool fixed)
{
    set<Node const*> ancestors;

#ifndef _RWSTD_NO_MEMBER_TEMPLATES
    ancestors.insert(gv->nodes().begin(), gv->nodes().end());
#else
    //Workaround for Solaris libCstd
    for (vector<StochasticNode const *>::const_iterator p = gv->nodes().begin();
         p != gv->nodes().end(); ++p)
      {
        ancestors.insert(*p);
      }
#endif 

    vector<DeterministicNode *> const &dnodes = gv->deterministicChildren();    
    for (unsigned int j = 0; j < dnodes.size(); ++j) {
	if (dnodes[j]->isClosed(ancestors, DNODE_POWER, fixed)) {
	    ancestors.insert(dnodes[j]);
	}
	else {
	    return false;
	}
    }
    
    return true;
}

} //namespace jags
