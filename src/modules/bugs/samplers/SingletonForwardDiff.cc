#include "SingletonForwardDiff.h"

#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <graph/DeterministicNode.h>

#include <map>

using std::map;
using std::vector;
using std::list;

namespace jags {

    static vector<unsigned long>
    child_indices(Node *node, map<Node const *, unsigned long> const &dmap) {
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
    
    namespace bugs {

	SingletonForwardDiff::SingletonForwardDiff(SingletonGraphView const *gv)
	    : _gv(gv)
	{
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
	    for (auto p = dn.begin(); p != dn.end(); ++p) {
		_edges2.push_back(child_indices(*p, dmap));
	    }

	    /*
	    //vector<vector<unsigned long>> _edges2;
	    _offsets.reserve(dn.size() + 1);
	    unsigned long L = 0;
	    _offsets.push_back(L);
	    for (auto p = dn.begin(); p != dn.end(); ++p) {
		L += (*p)->length();
		_offsets.push_back(L);
	    }
	    */
	    
	    /* Terminal nodes are the parents of the the stochastic children
	       in the GraphView. The sampled node may also be a terminal node.
	       In this case we use the index 0 and for deterministic nodes
	       we offset the index by 1 */
	    vector<StochasticNode *> const &sn = gv->stochasticChildren();
	    _nrow = 0;
	    for (auto p = sn.begin(); p != sn.end(); ++p) {
		vector<Node const *> const &parents = (*p)->parents();
		for (auto q = parents.begin(); q != parents.end(); ++q) {
		    if (*q == snode) {
			_terminal_nodes.push_back(0);
			_nrow += snode->length();
		    }
		    else {
			auto r = dmap.find(*q);
			if (r != dmap.end()) {
			    _terminal_nodes.push_back(r->second + 1);
			    _nrow += (*q)->length();
			}
		    }
		}
	    }
	}

	void SingletonForwardDiff::calCoef(double *coef, unsigned int chain) const
	{
	    StochasticNode *snode = _gv->node();
	    unsigned long nparam = snode->length();
	    vector<DeterministicNode *> const &dn = _gv->deterministicChildren();
	    
	    vector<vector<double>> grad(dn.size());
	    for (unsigned long i = 0; i < grad.size(); ++i) {
		grad[i] = vector<double>(nparam * dn[i]->length());
	    }

	    /*
	    for (unsigned long i = 0; i < _edges1.size(); ++i) {
		unsigned long j = _edges1[i];
		DeterministicNode const *c = dn[j];
	    }
	    for (auto r = _map.begin(); r != _map.end(); ++r) {
		vector<double> &v = r->second;
		fill(v.begin(), v.end(), 0.0);
	    }
	    */

	    // Edges from sampled node to its deterministic children
	    for (auto c = _edges1.begin(); c != _edges1.end(); ++c) {
		dn[*c]->gradient(grad[*c].data(), snode, chain);
	    }

	    // Edges between deterministic children
	    for (unsigned long p = 0; p < _edges2.size(); ++p) {
		for (auto c = _edges2[p].begin(); c != _edges2[p].end(); ++c) {
		
		    unsigned long nchild = dn[*c]->length();
		    unsigned long nparent = dn[p]->length();
		    vector<double> wk(nchild * nparent);
		    dn[*c]->gradient(wk.data(), dn[p], chain);

		    vector<double> &gc = grad[*c];
		    vector<double> &gp = grad[p];
		    for (unsigned long i = 0; i < nchild; ++i) {
			for (unsigned long j = 0; j < nparent; ++j) {
			    for (unsigned long k = 0; k < nparam; ++k) {
				gc[i + nchild*k] += wk[i + nchild*j] * gp[j + nparent*k];
			    }
			}
		    }
		}
	    }

	    /********
	    
	    list<DeterministicNode *> const *children = snode->deterministicChildren();
	    for (auto c = children->begin(); c != children->end(); ++c) {
		auto r = _map.find(*c);
		if (r != node_map.end()) {
		    (*c)->gradient(r->second, snode, chain);
		}
	    }
	    */

	    /*
	    // Edges between deterministic nodes
	    vector<DeterministicNode *> const &dn = _gv->deterministicChildren();
	    for (auto p = dn.begin(); p != dn.end(); ++p) {
		unsigned long n = (*p)->length();
		vector<Node *> const &children = (*p)->deterministicChildren();
		//vector<double> &v = _map.find(*p).second;
		for (auto c = children.begin(); c != children.end(); ++c) {
		    r = _map.find(*c);
		    if (r != node_map.end()) {
			vector<double> &v = r-second;
			unsigned long m = (*c)->length();
			vector<double> wk(m*n);
			(*c)->gradient(wk.data(), *p, chain);
			for (unsigned int i = 0; i < m; ++i) {
			    for (unsigned long j = 0; j < n; ++j) {
				for (unsigned long k = 0; k < n; ++k);
				r.second[i + n * j] += wk[i + m*k] * v[k + n*j]; //fixme
			    }
			}
		    }
		}
	    }
	    */

	    for (unsigned long i = 0; i < _terminal_nodes.size(); ++i) {
		unsigned long t = _terminal_nodes[i];
		if (t == 0) {
		    for (unsigned j = 0; j < nparam; ++j) {
			for (unsigned k = 0; k < nparam; ++k) {
			    coef[j + nparam*k] = (j==k) ? 1.0 : 0.0;
			}
		    }
		    coef += nparam * nparam;
		}
		else {
		    copy(grad[t-1].begin(), grad[t-1].end(), coef);
		    coef += dn[t-1]->length() * nparam;
		}
	    }

	    /*
		DeterministicNode *tnode = [_terminal_indices[i]];
		
	    }
	    
	    //Read off the gradients at the terminal nodes and concatenate them
	    vector<StochasticNode*> const &sc = _gv->stochasticChildren();
	    unsigned long L = 0;
	    for (auto p = sc.begin(); p != sc.end(); ++p) {
		for(auto q = (*p)->parents().begin(); q != (*p)->parents().end(); ++q) {
		    auto r = _map.find(*q);
		    if (r != _map.end()) {
			vector<double> const &v = r.second;
			copy(v.begin(), v.end(), grad + L);
			L += v.size();
		    }
		}
	    }

	    r = _map.find(_gv->node());
	    r.second[i] = 0;
	    */
	}

	unsigned long SingletonForwardDiff::nrow() const
	{
	    return _nrow;
	}

	unsigned long SingletonForwardDiff::ncol() const
	{
	    return _gv->length();
	}

    } /* Namespace bugs */
} /* Namespace jags */
