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
    
    namespace bugs {

	SingletonForwardDiff::SingletonForwardDiff(SingletonGraphView const *gv,
						   vector<Node const *> const &dependent_nodes)
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
	    //FIXME: Can simplify this by pushing more code into child_indicesx
	    _edges = vector<vector<unsigned long>>(dn.size());
	    vector<unsigned long> children = child_indices(snode, dmap);
	    for (auto c = children.begin(); c != children.end(); ++c) {
		_edges[*c].push_back(0);
	    }
	    for (unsigned int p = 0; p < dn.size(); ++p) {
		children = child_indices(dn[p], dmap);
		for (auto c = children.begin(); c != children.end(); ++c) {
		    _edges[*c].push_back(p+1);
		}
	    }
	    
	    _nrow = 0;
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
		}
	    }

	    _workspace_nrow = 0;
	    for (auto p = dn.begin(); p != dn.end(); ++p) {
		_workspace_nrow += (*p)->length();
	    }
	}

	void SingletonForwardDiff::calCoef(double *coef, unsigned int chain, bool transpose) const
	{
	    StochasticNode *snode = _gv->node();
	    unsigned long nparam = snode->length();
	    vector<DeterministicNode *> const &dn = _gv->deterministicChildren();

	    // A common workspace hods the gradient matrices for each deterministic node
	    vector<double> gradient_workspace(_workspace_nrow * nparam);
	    /* Vector of pointers to the start of the gradient array for each node.
	       The pointers are constructed on the fly */
	    vector<double *> gradient_start(dn.size() + 1);
	    gradient_start[0] = gradient_workspace.data();
	       
	    // Calculate gradient for each edge
	    vector<double> wk; // Reusable workspace for edges
	    for (unsigned long c = 0; c < dn.size(); ++c) {

		unsigned long nchild = dn[c]->length();
		double *gc = gradient_start[c]; // nchild x nparam

		for (auto r = _edges[c].begin(); r != _edges[c].end(); ++r) {

		    if (*r == 0) {
			// Parent is sampled node
			dn[c]->gradient(gc, snode, chain);
		    }
		    else {
			// Parent is deterministic node
			unsigned long p = *r - 1;
			unsigned long nparent = dn[p]->length();
			// Reset workspace
			wk.resize(nchild * nparent);
			fill(wk.begin(), wk.end(), 0);
			// Get the gradient for the edge p -> c
			dn[c]->gradient(wk.data(), dn[p], chain);
			// Apply chain rule to get gradient wrt sampled node
			double *gp = gradient_start[p]; // nparent x nparam
			//FIXME Replace this with a call to dgemm
			for (unsigned long i = 0; i < nchild; ++i) {
			    for (unsigned long j = 0; j < nparam; ++j) {
				for (unsigned long k = 0; k < nparent; ++k) {
				    gc[i + nchild*j] += wk[i + nchild*k] * gp[k + nparent*j];
				}
			    }
			}
		    }
		}
		//Define start of gradient matrix for dn[c+1]
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
