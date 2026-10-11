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
	    _edges1 = child_indices(snode, dmap);
	    _edges2 = vector<vector<unsigned long>>(dn.size());
	    for (unsigned int p = 0; p < dn.size(); ++p) {
		vector<unsigned long> children = child_indices(dn[p], dmap);
		for (auto c = children.begin(); c != children.end(); ++c) {
		    _edges2[*c].push_back(p);
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

	    unsigned long gnrow = 0;
	    _gindex.reserve(dn.size() + 1);
	    _gindex.push_back(0);
	    for (auto p = dn.begin(); p != dn.end(); ++p) {
		gnrow += (*p)->length();
		_gindex.push_back(gnrow);
	    }
	}

	void SingletonForwardDiff::calCoef(double *coef, unsigned int chain, bool transpose) const
	{
	    StochasticNode *snode = _gv->node();
	    unsigned long nparam = snode->length();
	    vector<DeterministicNode *> const &dn = _gv->deterministicChildren();

	    // A common workspace hods the gradient matrices for each deterministic node
	    vector<double> gradient_workspace(_gindex.back() * nparam);
	    auto get_gradmat = [&](unsigned long i) {
		// Returns a pointer to the start of the gradient
		// array for the deterministic node with index i
		return gradient_workspace.data() + _gindex[i] * nparam;
	    };
	    
	    // Edges from sampled node to deterministic children
	    for (auto c = _edges1.begin(); c != _edges1.end(); ++c) {
		dn[*c]->gradient(get_gradmat(*c), snode, chain);
	    }

	    // Edges between deterministic children
	    vector<double> wk; // Reusable workspace for edges
	    for (unsigned long c = 0; c < dn.size(); ++c) {

		unsigned long nchild = dn[c]->length();
		double *gc = get_gradmat(c); // nchild x nparam

		for (auto p = _edges2[c].begin(); p != _edges2[c].end(); ++p) {

		    unsigned long nparent = dn[*p]->length();
		    // Reset workspace
		    wk.resize(nchild * nparent);
		    fill(wk.begin(), wk.end(), 0);
		    // Get the gradient for the edge p -> c
		    dn[c]->gradient(wk.data(), dn[*p], chain);
		    // Apply chain rule to get gradient wrt sampled node
		    double *gp = get_gradmat(*p); // nparent x nparam
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
		    double *gt = get_gradmat(t-1);
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
