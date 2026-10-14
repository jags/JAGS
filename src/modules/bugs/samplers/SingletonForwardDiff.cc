#include "SingletonForwardDiff.h"

#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <graph/DeterministicNode.h>

#include <map>

using std::map;
using std::vector;
using std::list;

//debuggin
#include <module/ModuleError.h>

namespace jags {

    /* Utility function for building adjacency matrix */
    static void
    push_child_indices(Node *node, map<Node const *, unsigned long> const &dmap,
		       vector<vector<unsigned long>> &edges, unsigned long i) {
	list<DeterministicNode *> const *children = node->deterministicChildren();
	for (auto p = children->begin(); p != children->end(); ++p) {
	    auto r = dmap.find(*p);
	    if (r != dmap.end()) {
		edges[r->second].push_back(i);
	    }
	}
    }

    unsigned long length(Node const *node) { return node->length(); }

    template<typename T> unsigned long length(vector<T> const &v) { return v.size(); }
    
    template<typename T> unsigned long sumLengths(vector<T> const &v) {
	unsigned long L = 0;
	for (auto p = v.begin(); p != v.end(); ++p) {
	    L += length(*p);
	}
	return L;
    }
    
    namespace bugs {

	SingletonForwardDiff::SingletonForwardDiff(SingletonGraphView const *gv,
						   vector<Node const *> const &dependent_nodes)
	    : _gv(gv), _nrow(0), _workspace_nrow(sumLengths(gv->deterministicChildren()))
	{
	    StochasticNode *snode = gv->node();
	    vector<DeterministicNode *> const &dn = gv->deterministicChildren();

	    //debuggin sanity check
	    for (auto p = dependent_nodes.begin(); p != dependent_nodes.end(); ++p) {
		if ((*p) == nullptr)
		    throwLogicError("Oh no!");
	    }
	    
	    /* We need a temporary map to replace pointers to the deterministic
	       children with their indices in the vector dn. */
	    map<Node const *, unsigned long> dmap;
	    for (unsigned long i = 0; i < dn.size(); ++i) {
		dmap.insert({dn[i], i});
	    }

	    /* Create adjacency matrix of edges between deterministic
	       nodes, first as nested vectors ... */
	    vector<vector<unsigned long>> edges(dn.size());
	    push_child_indices(snode, dmap, edges, 0);
	    for (unsigned int p = 0; p < dn.size(); ++p) {
		push_child_indices(dn[p], dmap, edges, p+1);
	    }
	    /* ... then convert to compressed column format */
	    _edges_p.reserve(dn.size() + 1); // column indices
	    _edges_i.reserve(sumLengths(edges)); // row numbers
	    _edges_p.push_back(0);
	    for (auto p = edges.begin(); p != edges.end(); ++p) {
		_edges_p.push_back(_edges_p.back() + p->size());
		_edges_i.insert(_edges_i.end(), p->begin(), p->end());
	    }

	    /* Convert dependent nodes to indices, dropping nodes with trivial zero gradient  */
	    for (auto p = dependent_nodes.begin(); p != dependent_nodes.end(); ++p) {
		if (*p == snode) {
		    _dependent_nodes.push_back(0);
		    _nrow += snode->length();
		}
		else {
		    auto r = dmap.find(*p);
		    if (r != dmap.end()) {
			_dependent_nodes.push_back(r->second + 1);
			_nrow += (*p)->length();
		    }
		    // Else dropped due to zero gradient
		}
	    }
	}

	void SingletonForwardDiff::calCoef(double *coef, unsigned int chain, bool transpose) const
	{
	    StochasticNode *snode = _gv->node();
	    unsigned long nparam = snode->length();
	    vector<DeterministicNode *> const &dn = _gv->deterministicChildren();

	    // A common workspace hods the gradient matrices for each deterministic node
	    vector<double> dnode_workspace(_workspace_nrow * nparam);
	    // Vector of pointers to the start of the gradient array for each deterministic node.
	    vector<double *> gradient_start(dn.size() + 1);
	    gradient_start[0] = dnode_workspace.data();
	       
	    // Calculate gradient for each edge
	    vector<double> edge_workspace; // Reusable workspace for edges
	    for (unsigned long c = 0; c < dn.size(); ++c) {

		unsigned long nchild = dn[c]->length();
		double *gc = gradient_start[c]; // nchild x nparam

		for (auto r = _edges_p[c]; r < _edges_p[c+1]; ++r) {
		    if (_edges_i[r] == 0) {
			// Parent is sampled node
			dn[c]->gradient(gc, snode, chain);
		    }
		    else {
			// Parent is deterministic node
			unsigned long p = _edges_i[r] - 1;
			unsigned long nparent = dn[p]->length();
			// Reset workspace
			edge_workspace.resize(nchild * nparent); // nchild x nparent
			fill(edge_workspace.begin(), edge_workspace.end(), 0);
			// Get the gradient for the edge p -> c
			dn[c]->gradient(edge_workspace.data(), dn[p], chain);
			// Apply chain rule to get gradient wrt sampled node
			double *gp = gradient_start[p]; // nparent x nparam
			//FIXME Replace this with a call to dgemm
			for (unsigned long i = 0; i < nchild; ++i) {
			    for (unsigned long j = 0; j < nparam; ++j) {
				for (unsigned long k = 0; k < nparent; ++k) {
				    gc[i + nchild*j] += edge_workspace[i + nchild*k] * gp[k + nparent*j];
				}
			    }
			}
		    }
		}
		//Define start of gradient matrix for next deterministic child dn[c+1]
		gradient_start[c+1] = gradient_start[c] + nchild * nparam;
	    }
    
	    /* Stack gradient matrices for dependent nodes to form coefficient matrix
	       Some attention is required here as the gradient matrices cannot
	       simply be concatenated.
	    */
	    unsigned long j0 = 0;
	    for (unsigned long i = 0; i < _dependent_nodes.size(); ++i) {
		unsigned long t = _dependent_nodes[i];
		if (t == 0) {
		    //Dependent node is the sampled node
		    for (unsigned long j = 0; j < nparam; ++j) {
			for (unsigned long k = 0; k < nparam; ++k) {
			    coef[j0 + j + _nrow*k] = (j==k) ? 1.0 : 0.0;
			}
		    }
		    j0 += nparam;
		}
		else {
		    //Dependent node is deterministic node dn[t-1]
		    unsigned long L = dn[t-1]->length();
		    double *gt = gradient_start[t-1];
		    for (unsigned long j = 0; j < L; ++j) {
			for (unsigned long k = 0; k < nparam; ++k) {
			    coef[j0 + j + _nrow*k] = gt[j + L*k];
			}
		    }
		    j0 += L;
		}
	    }
	    
	    if (transpose) {
		vector<double> tcoef(_nrow * nparam);
		unsigned long nparam = _gv->length();
		for (unsigned long j = 0; j < nparam; ++j) {
		    for (unsigned long k = 0; k < _nrow; ++k) {
			tcoef[j + nparam*k] = coef[k + _nrow *j];
		    }
		}
		copy(tcoef.begin(), tcoef.end(), coef);
	    }

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
