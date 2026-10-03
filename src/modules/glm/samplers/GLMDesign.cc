#include "GLMDesign.h"

#include <sampler/SingletonGraphView.h>
#include <sampler/Linear.h>
#include <graph/StochasticNode.h>
#include <graph/DeterministicNode.h>
#include <graph/LinkNode.h>

#include <map>
#include <stdexcept>

using std::map;
using std::vector;
using std::list;
using std::invalid_argument;

//#include <algorithm>

namespace jags {

    static vector<unsigned long>
    child_indices(Node *node, map<Node const *, unsigned long> const &dmap) {
	/* Convert deterministic children of node to a vector of
	 * numeric indices using the supplied map */
	vector<unsigned long> ans;
	list<DeterministicNode *> const *children = node->deterministicChildren();
	for (auto p = children->begin(); p != children->end(); ++p) {
	    auto r = dmap.find(*p);
	    if (r != dmap.end()) {
		ans.push_back(r->second);
	    }
	}
	return ans;
    }
    
    namespace glm {

	GLMDesign::GLMDesign(SingletonGraphView const *gv,
				       vector<StochasticNode *> const &stoch_children)
	    : _gv(gv)
	{
	    if (!checkLinear(gv, false, true)) {
		throw invalid_argument("GraphView is not a GLM");
	    }
	    
	    // Access data from gv
	    StochasticNode *snode = gv->node();
	    vector<DeterministicNode *> const &dn = gv->deterministicChildren();
	    
	    /* We need a temporary map to replace pointers to the deterministic
	       children with their indices in the vector dn. */
	    map<Node const *, unsigned long> dmap;
	    for (unsigned long i = 0; i < dn.size(); ++i) {
		dmap.insert({dn[i], i});
	    }

	    /* There are two kinds of edge. The first is from the sampled node
	       (snode) into the deterministic children(dn). The second is between
	       deterministic children. We keep separate records of these edges
	       as they are handled differently in calCoef */
	    _edges1 = child_indices(snode, dmap);
	    _edges2 = vector<vector<unsigned long>>(dn.size()); //FIXME should be in initialization
	    for (unsigned int p = 0; p < dn.size(); ++p) {
		vector<unsigned long> ichildren = child_indices(dn[p], dmap);
		for (auto c = ichildren.begin(); c != ichildren.end(); ++c) {
		    _edges2[*c].push_back(p);
		}
	    }
	    
	    /* Terminal nodes are the parents of the the stochastic
	       children supplied to the constructor.  As this is a GLM
	       we take only the first parent. If this is a link node
	       then we take the grandparent node to get the linear
	       predictor.
	       
	       The sampled node may also be a terminal node.  In this
	       case we use the index 0 and for deterministic nodes we
	       offset the index by 1
	    */
	    _nrow = 0;
	    for (auto p = stoch_children.begin(); p != stoch_children.end(); ++p) {
		Node const *tnode = (*p)->parents()[0];
		if (dynamic_cast<LinkNode const *>(tnode) != NULL) {
		    tnode = tnode->parents()[0];
		}
		if (tnode == snode) {
		    _terminal_nodes.push_back(0);
		    _nrow += snode->length();
		}
		else {
		    auto r = dmap.find(tnode);
		    if (r != dmap.end()) {
			_terminal_nodes.push_back(r->second + 1);
			_nrow += tnode->length();
		    }
		}
	    }
	}

	void GLMDesign::calCoef(double *coef, unsigned int chain) const
	{
	    StochasticNode *snode = _gv->node();
	    unsigned long nparam = snode->length();
	    vector<DeterministicNode *> const &dn = _gv->deterministicChildren();
	    
	    vector<vector<double>> grad(dn.size());
	    for (unsigned long i = 0; i < grad.size(); ++i) {
		grad[i] = vector<double>(nparam * dn[i]->length());
	    }

	    // Edges from sampled node to its deterministic children
	    for (auto c = _edges1.begin(); c != _edges1.end(); ++c) {
		dn[*c]->gradient(grad[*c].data(), snode, chain);
	    }

	    // Edges between deterministic children
	    for (unsigned long c = 0; c < _edges2.size(); ++c) {
		for (auto p = _edges2[c].begin(); p != _edges2[c].end(); ++p) {
		
		    unsigned long nchild = dn[c]->length();
		    unsigned long nparent = dn[*p]->length();
		    vector<double> wk(nchild * nparent); // nchild x nparent
		    dn[c]->gradient(wk.data(), dn[*p], chain);

		    //FIXME Replace this with a call to dgemm
		    vector<double> &gc = grad[c]; // nchild x nparam
		    vector<double> &gp = grad[*p]; // nparent x nparam
		    for (unsigned long i = 0; i < nchild; ++i) {
			for (unsigned long j = 0; j < nparam; ++j) {
			    for (unsigned long k = 0; k < nparent; ++k) {
				gc[i + nchild*j] += wk[i + nchild*k] * gp[k + nparent*j];
				// nchild x nparam   nchild x nparent   nparent x nparam
			    }
			}
		    }
		}
	    }

	    /* Stack gradient matrices for terminal nodes form coefficient matrix
	       Some attention is required here as the gradient matrices cannot
	       simply be concatenated.
	     */
	    unsigned long j0 = 0;
	    for (unsigned long i = 0; i < _terminal_nodes.size(); ++i) {
		unsigned long t = _terminal_nodes[i];
		if (t == 0) {
		    //Terminal node is the sampled node
		    for (unsigned long j = 0; j < nparam; ++j) {
			for (unsigned long k = 0; k < nparam; ++k) {
			    coef[j0 + j + _nrow*k] = (j==k) ? 1.0 : 0.0;
			}
		    }
		    j0 += nparam;
		}
		else {
		    //Terminal node is deterministic node dn[t-1]
		    unsigned long L = dn[t-1]->length();
		    for (unsigned long j = 0; j < L; ++j) {
			for (unsigned long k = 0; k < nparam; ++k) {
			    coef[j0 + j + _nrow*k] = grad[t-1][j + L*k];
			}
		    }
		    j0 += L;
		}
	    }
	}

	unsigned long GLMDesign::nrow() const
	{
	    return _nrow;
	}

	unsigned long GLMDesign::ncol() const
	{
	    return _gv->length();
	}
	
    } /* Namespace glm */
} /* Namespace jags */
