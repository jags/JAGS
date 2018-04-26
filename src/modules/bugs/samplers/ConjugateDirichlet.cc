#include <config.h>

#include "ConjugateDirichlet.h"

#include <rng/RNG.h>
#include <graph/MixtureNode.h>
#include <graph/StochasticNode.h>
#include <graph/AggNode.h>
#include <graph/NodeError.h>
#include <sampler/SingletonGraphView.h>
#include <module/ModuleError.h>

#include <set>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <climits>

#include <JRmath.h>

using std::vector;
using std::set;
using std::string;
using std::copy;
using std::map;
using std::pair;
using std::list;

namespace jags {

    Node const *findUniqueParent(Node const *node,
				 set<Node const *> const &nodeset)
    {
	/*
	  Utility function called by checkAggNode and checkMixNode.
	  
	  If node has a single parent in nodeset, then a pointer to
	  the parent is returned.  If there are multiple parents in
	  nodeset then a NULL pointer is returned.
	      
	  If no parents are in nodeset, a logic error is thrown.
	*/
	vector<Node const *> const &par = node->parents();
	Node const *param = 0;
	    
	for (unsigned int i = 0; i < par.size(); ++i) {
	    if (nodeset.count(par[i])) {
		if (param) {
		    if (param != par[i]) return 0;
		}
		else {
		    param = par[i];
		}
	    }
	}
	if (param == 0) {
	    throwLogicError("Error in ConjugateDirichlet::canSample");
	}
	    
	return param;
    }


    bool checkAggNode(AggNode const *anode, 
		      set<Node const *> const &nodeset)
    {
	/* 
	   Utility function called by ConjugateDirichlet::canSample
	       
	   The aggregate node (anode) must have only one parent in
	   nodeset.  This parent must be embedded entirely within
	   anode (i.e. we cannot take a subset of it) and the order of
	   the elements cannot be permuted.
	*/
	    
	Node const *param = findUniqueParent(anode, nodeset);
	if (param == 0) return false;
	    
	//Check that parent is entirely contained in anode with
	//offsets in ascending order

	vector<Node const *> const &par = anode->parents();
	vector<unsigned long> const &off = anode->offsets();

	unsigned int j = 0;
	for (unsigned int i = 0; i < par.size(); ++i) {
	    if (par[i] == param && off[i] != j++) return false;
	}
	if (j != param->length()) return false;
	    
	return true;
    }

    bool checkMixNode(MixtureNode const *mnode, 
		      set<Node const *> const &nodeset)
    {
	/*
	  Utility function called by ConjugateDirichlet::canSample
	      
	  None of the indices may depend on nodeset. Among the other
	  parents, only 1 may be in nodeset (although it may appear
	  several times)
	*/
	    
	vector<Node const *> const &par = mnode->parents();
	    
	//Check indices
	unsigned long nindex = mnode->index_size();
	for (unsigned long i = 0; i < nindex; ++i) {
	    if (nodeset.count(par[i])) return false;
	}
	    
	//Check for unique parent
	return findUniqueParent(mnode, nodeset) != 0;
    }

    bool isMix(SingletonGraphView const *gv)
    {
	/* 
	   Utility function called by constructor. It returns true if
	   any of the deterministic descendants are mixture nodes.
	   
	   When there are no mixture nodes, the update calculations
	   can be simplified.
	   
	   This static function is called by the constructor and the
	   result is stored in the variable _mix.
	*/
	vector<DeterministicNode*> const &dchild = gv->deterministicChildren();
	for (unsigned int i = 0; i < dchild.size(); ++i) {
	    if (isMixture(dchild[i])) return true;
	}
	return false;
    }

    vector<unsigned long> makeTree(SingletonGraphView const *gv)
    {
	/* 
	   If canSample is true then the nodes in the GraphView gv
	   form a tree with the sampled node as its root.
	       
	   This function creates an integer vector "tree" denoting the
	   parents of the deterministic descendants: "tree[i] = j" for
	   0 < j < i means that the unique parent of deterministic
	   node "i" in the GraphView is deterministic node "j".  If
	   the unique parent is the sampled node then we denote this
	   by "tree[i] = ULONG_MAX".

	   This static function is called by the constructor, and the
	   result is stored in the variable _tree;
	*/
	vector<DeterministicNode*> const &dchild = gv->deterministicChildren();
	StochasticNode *snode = gv->node();

	vector<unsigned long> tree(dchild.size(), ULONG_MAX);
	    
	set<Node const *> nodeset;
	nodeset.insert(snode);
	    
	for (unsigned int j = 0; j < dchild.size(); ++j) {
	    Node const *parent = findUniqueParent(dchild[j], nodeset);
	    if (parent == 0) {
		throwLogicError("Invalid tree in ConjugateDirichlet");
	    }
	    if (parent != snode) {
		for (unsigned int k = 0; k < j; ++k) {
		    if (parent == dchild[k]) {
			tree[j] = k;
			break;
		    }
		}
		if (tree[j] == ULONG_MAX) {
		    throwLogicError("Invalid tree in ConjugateDirichlet");
		}
	    }
	    nodeset.insert(dchild[j]);
	}	
	    
	return tree;
    }

    vector<vector<unsigned long> > makeOffsets(SingletonGraphView const *gv,
				       vector<unsigned long> const &tree)
    {
	vector<DeterministicNode *> const &dchild = gv->deterministicChildren();
	vector<vector<unsigned long> > offsets(dchild.size());
	StochasticNode const *snode = gv->node();

	for (unsigned int i = 0; i < dchild.size(); ++i) {
	
	    unsigned long j = tree[i]; //index of parent in dchild
	    
	    if (isMixture(dchild[i])) {
		if (j != ULONG_MAX) offsets[i] = offsets[j];
	    }
	    else if (AggNode const *a = dynamic_cast<AggNode const*>(dchild[i]))
	    {
		vector<Node const *> const &par_i = a->parents();
		vector<unsigned long> const &off_i = a->offsets();
		
		Node const *target = (j == ULONG_MAX) ? 
		    static_cast<Node const *>(snode) : 
		    static_cast<Node const *>(dchild[j]);
		
		if (j == ULONG_MAX || offsets[j].empty()) {
		    for (unsigned int k = 0; k < par_i.size(); ++k) {
			if (par_i[k] == target) {
			    offsets[i].push_back(k);
			}
		    }
		}
		else {
		    unsigned int p = 0;
		    for (unsigned int k = 0; k < par_i.size(); ++k) {
			if (par_i[k] == target && off_i[k] == offsets[j][p]) {
			    offsets[i].push_back(k);
			    p++;
			}
		    }
		}
		if (offsets[i].size() != snode->length()) {
		    throwLogicError("Offset error in ConjugateDirichlet");
		}
	    }
	    else {
		throwLogicError("Invalid child in ConjugateDirichlet");
	    }
	}
	return offsets;
    }
	
    namespace bugs {
	
	bool 
	ConjugateDirichlet::canSample(StochasticNode *snode, Graph const &graph)
	{
	    /*
	      The ConjugateDirichlet sampler aims to be as general as
	      possible, while minimizing the calculations that have to
	      be done at each iteration.
      
	      The stochastic children may be multinomial or
	      categorical, and the deterministic children may include
	      mixture nodes and aggregate nodes. However there are
	      some restrictions that are implemented in the
	      checkAggNode and checkMixNode functions.
      
	      1) The deterministic descendents must form a tree, so
	      that there is only one pathway from the sample node
	      (snode) to each stochastic child.
      
	      2) Aggregate nodes can only be used to embed snode (or
	      one of its determistic descendants) in a larger node,
	      not to take a subset.

	      3) When a node X is embedded in an aggregate node Y, the
	      order of the elements of X cannot be permuted in Y.

	      These restrictions ensure that the _off values used by
	      the update method are the same across chains and across
	      iterations, and can therefore be calculated once in the
	      constructor.

	      Fascinating fact: Neither dmulti nor dcat requires a
	      normalized probability argument.  This allows us to
	      embed snode inside a larger aggregate node without
	      worrying about what the other values are. For example:
	      
	      Y ~ dcat(p)
	      p[1:3] ~ ddirch(alpha)
	      p[4:5] <- W
	      
	      We can completely ignore the values of W (p[4] and
	      p[5]). Normally these will be zero but we don't require
	      it.
	    */
    
	    if (getDist(snode) != DIRCH)
		return false;
	    
	    if (isBounded(snode))
		return false;
    
	    SingletonGraphView gv(snode, graph);
	    vector<DeterministicNode*> const &dchild = gv.deterministicChildren();
	    vector<StochasticNode *> const &schild = gv.stochasticChildren();
    
	    // Check stochastic children
	    for (unsigned int i = 0; i < schild.size(); ++i) {
		vector<Node const *> const &param = schild[i]->parents();
		if (isBounded(schild[i])) {
		    return false; //Truncated
		}
		switch(getDist(schild[i])) {
		case CAT:
		    break;
		case MULTI:
		    if (gv.isDependent(param[1])) return false;
		    break;
		default:
		    return false;
		}
	    }

	    // Check deterministic descendants
	    set<Node const *> nodeset;
	    nodeset.insert(snode);
	    for (unsigned int j = 0; j < dchild.size(); ++j) {
		if (MixtureNode const *m = dynamic_cast<MixtureNode const*>(dchild[j]))
		{
		    if (!checkMixNode(m, nodeset)) return false;	    
		}
		else if (AggNode const *a = dynamic_cast<AggNode const*>(dchild[j]))
		{
		    if (!checkAggNode(a, nodeset)) return false;
		}
		else {
		    return false;
		}
		nodeset.insert(dchild[j]);
	    }
	    
	    return true;
	}
	
	ConjugateDirichlet::ConjugateDirichlet(SingletonGraphView const *gv)
	    : ConjugateMethod(gv), _mix(isMix(gv)), _tree(makeTree(gv)),
	      _offsets(gv->stochasticChildren().size()),
	      _leaves(gv->stochasticChildren().size(), ULONG_MAX)
	{
	    vector<vector<unsigned long> > offsets = makeOffsets(gv, _tree);
	    
	    //Create a map from stochastic children onto their indices
	    vector<StochasticNode *> const &schild =
		gv->stochasticChildren();

	    map<StochasticNode const *, unsigned long> smap;
	    for (unsigned long s = 0; s < schild.size(); ++s) {
		smap.insert(pair<StochasticNode const*, unsigned long>(schild[s], s));
	    }

	    vector<DeterministicNode *> const &dchild = 
		gv->deterministicChildren();
	    for (unsigned long d = 0; d < dchild.size(); ++d) {
	
		list<StochasticNode *> const *dc = 
		    dchild[d]->stochasticChildren();

		for (list<StochasticNode *>::const_iterator q = dc->begin(); 
		     q != dc->end(); ++q) 
		{
		    map<StochasticNode const *, unsigned long >::iterator r = smap.find(*q);
		    if (r != smap.end()) {
			_leaves[r->second] = d;
			_offsets[r->second] = offsets[d];
			smap.erase(r);
		    }
		}

	    }
	}

	bool ConjugateDirichlet::isActiveLeaf(unsigned long i,
					      unsigned int chain) const
	{
	    if (!_mix) return true;
	    return isActiveTree(_leaves[i], chain);
	}
	
	bool ConjugateDirichlet::isActiveTree(unsigned long i,
					      unsigned int chain) const 
	{
	    /*
	      Returns true if there is an active path from the sampled
	      node to the deterministic node with index "i".  A path
	      is not active if it is blocked by a mixture node that is
	      set to another mixture component.
	    */

	    if (i == ULONG_MAX) {
		return true; //We have reached the sampled node
	    }

	    vector<DeterministicNode*> const &dchild = 
		_gv->deterministicChildren();
    
	    if (MixtureNode const *m = asMixture(dchild[i])) {
		if (_tree[i] == ULONG_MAX) {
		    if (m->activeParent(chain) != _gv->node())
			return false;
		}
		else {
		    if (m->activeParent(chain) != dchild[_tree[i]])
			return false;
		}
	    }
    
	    return isActiveTree(_tree[i], chain);
	}

void ConjugateDirichlet::update(unsigned int chain, RNG *rng) const
{
    StochasticNode *snode = _gv->node();
    unsigned long size = snode->length();
    double *alpha = new double[size];
    double *xnew = new double[size];

    double const *prior = snode->parents()[0]->value(chain);
    for (unsigned long i = 0; i < size; ++i) {
	alpha[i] = prior[i];
    }

    vector<StochasticNode *> const &schild = _gv->stochasticChildren();
    for (unsigned int i = 0; i < schild.size(); ++i) {
	unsigned int index = 0;
	double const *N = 0;

	if (isActiveLeaf(i, chain)) {
	    switch(_child_dist[i]) {
	    case MULTI:
		N = schild[i]->value(chain);
		if (_offsets[i].empty()) {
		    for (unsigned int j = 0; j < size; ++j) {
			alpha[j] += N[j];
		    }
		}
		else {
		    for (unsigned int j = 0; j < size; ++j) {
			alpha[j] += N[_offsets[i][j]];
		    }
		}
		break;
	    case CAT:
		index = static_cast<unsigned int>(*schild[i]->value(chain)) - 1;
		if (_offsets[i].empty()) {
		    alpha[index] += 1;
		}
		else {
		    for (unsigned int j = 0; j < size; ++j) {
			if (index == _offsets[i][j]) {
			    alpha[j] += 1;
			    break;
			}
		    }
		}
		break;
	    default:
		throwLogicError("Invalid distribution in ConjugateDirichlet");
	    }
	}
    }

    /* Check structural zeros */
    for (unsigned int i = 0; i < size; ++i) {
	if (prior[i] == 0 && alpha[i] != 0) {
	    throwNodeError(snode, "Invalid likelihood for Dirichlet distribution with structural zeros");
	}
    }
  
    /* 
       Draw Dirichlet sample by drawing independent gamma random
       variates and then normalizing
    */

    double xsum = 0.0;
    for (unsigned long i = 0; i < size; ++i) {
	if (alpha[i] > 0) {
	    xnew[i] = rgamma(alpha[i], 1, rng);
	    xsum += xnew[i];
	}
	else {
	    xnew[i] = 0;
	}
    }
    for (unsigned long i = 0; i < size; ++i) {
	xnew[i] /= xsum;
    }

    if (_mix) {

	//Override default updating of deterministic nodes in gv

	snode->setValue(xnew, size, chain);

	vector<DeterministicNode*> const &dchild = _gv->deterministicChildren();
	vector<bool> modified(dchild.size(), false);
	for (unsigned long d = 0; d < dchild.size(); ++d) {
	    if (_tree[d] == ULONG_MAX) {
		MixtureNode *m = dynamic_cast<MixtureNode*>(dchild[d]);
		if (m == 0 || m->activeParent(chain) == snode) {
		    dchild[d]->deterministicSample(chain);
		    modified[d] = true;
		}
	    }
	    else if (modified[_tree[d]]) {
		MixtureNode *m = dynamic_cast<MixtureNode*>(dchild[d]);
		if (m == 0 || m->activeParent(chain) == dchild[_tree[d]]) {
		    dchild[d]->deterministicSample(chain);
		    modified[d] = true;
		}
	    }
	}
    }
    else {
	_gv->setValue(xnew, size, chain);
    }

    delete [] xnew;
    delete [] alpha;
}


    } //namespace jags
} //namespace bugs

