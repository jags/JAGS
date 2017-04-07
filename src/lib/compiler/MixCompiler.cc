/* Compiler rules for Mixture nodes.

   Node * getMixtureNode(ParseTree const * var, Compiler *compiler);

   This is painfully complicated, but the aim is to be as efficient as
   possible in the allocation of Mixture Nodes.

   Suppose we have an expression X[Y], where X is a vector node of length
   N, and Y is a discrete-valued scalar node. We know that Y can only
   legally take values 1, ..., N. So the simplest approach is to create
   a Mixture Node that takes value X[i] when Y = i, for i = 1, ... N.
   This mixture node has N+1 parents: the index node Y and the N subset
   nodes X[1], ... X[N].

   It may not be efficient to do this in cases where the range of Y is
   restricted.  For example, suppose Y only takes values in the range M1,
   ..., M2, where M1 > 1 and M2 < N.  Then the parents X[1], ... X[M1 -
   1] and X[M2 + 1] ... X[N] are redundant. This may be harmless, but
   there are circumstances under which it will lead to a compiler error.

   There are limited circumstances under which we can calculate the
   possible values of the index node. Currently, this is restricted to
   the case where all the Stochastic parents of Y are discrete
   and univariate. In this case, the function getMixtureNode1 is used
   to create an efficient Mixture Node. If this fails, we fall back on
   the default getMixtureNode2 function.
*/

#include <compiler/Compiler.h>
#include <compiler/ParseTree.h>
#include <model/SymTab.h>
#include <graph/StochasticNode.h>
#include <model/NodeArray.h>
#include <graph/MixtureNode.h>
#include <graph/GraphMarks.h>
#include <graph/NodeError.h>
#include <sampler/Sampler.h>
#include <util/nainf.h>

#include "MixCompiler.h"

#include <stdexcept>
#include <string>
#include <set>
#include <climits>
#include <algorithm>
#include <list>
#include <map>

using std::vector;
using std::pair;
using std::logic_error;
using std::runtime_error;
using std::string;
using std::set;
using std::map;
using std::copy;
using std::reverse;
using std::list;

namespace jags {

    //Saves the values of the nodes in chain 0 to the values vector
    //(presumed empty when the function is called)
    template<class T>
    void save(vector<T*> const &nodes, vector<vector<double> > &values)
    {
	for (typename vector<T*>::const_iterator p = nodes.begin();
	     p != nodes.end(); ++p)
	{
	    double const *pv = (*p)->value(0);
	    unsigned int n = (*p)->length();
	    vector<double> v(n);
	    copy(pv, pv + n, v.begin());
	    values.push_back(v);
	}
    }

    //Restores saved values to the nodes in chain 0
    template<class T>
    void restore(vector<T*> const &nodes, vector<vector<double> > const &values)
    {
	for (unsigned int j = 0; j < nodes.size(); ++j) {
	    vector<double> const &pv = values[j];
	    nodes[j]->setValue(&pv[0], pv.size(), 0);
	}
    }

//Structure to hold subset indices
struct SSI {
  Node const *node;
  vector<int> indices;
  SSI() : node(0) {};
};


/* 
   Given a vector of subset indices, this function modifies the argument
   "subsets", so that it contains a vector of pairs:
   - Index: Possible value of the variable indices
   - Range: Corresponding range of the array to take as a subset

   All possible values of Index are included in the vector, based on the
   the range of the array we are taking subsets from (default_range).
*/
static void getSubsetRanges(vector<pair<vector<int>, Range> > &subsets,  
			    vector<SSI> const &limits,
			    SimpleRange const &default_range)
{
    unsigned int ndim = limits.size();

    // Create upper and lower bounds
    vector<int> var_offset, var_lower, var_upper;
    vector<vector<int> > scope(ndim);
    for (unsigned int j = 0; j < ndim; ++j) {
	if (limits[j].node) {
	    var_offset.push_back(j);
	    var_lower.push_back(default_range.lower()[j]);
	    var_upper.push_back(default_range.upper()[j]);
	}
	else {
	    scope[j] = limits[j].indices;
	}
    }

    SimpleRange var_range(var_lower, var_upper); //range of variable indices
    for (RangeIterator p(var_range); !p.atEnd(); p.nextLeft()) {
	for (unsigned int k = 0; k < var_offset.size(); ++k) {
	    scope[var_offset[k]] = vector<int>(1, p[k]);
	}
	subsets.push_back(pair<vector<int>, Range>(p, Range(scope)));
    }
}

/*  Recursive function that returns true if node, or one of its
    deterministic descendants, is in target_set */
static bool hasDescendant(DeterministicNode *node, 
			  set<Node const*> const &target_set,
			  vector<DeterministicNode*> &dnodes,
			  map<DeterministicNode*, bool> &known_dnodes)
{
    //Skip previously visited nodes
    map<DeterministicNode*,bool>::iterator i = known_dnodes.find(node);
    if (i != known_dnodes.end()) {
	return i->second;
    }

    //The current node might be in the target set.
    bool ans = target_set.count(node);

    //Search all deterministic descendants. Even if node is in
    //target_set we still need to visit the descendants to ensure that
    //nodes are pushed onto the vector dnodes *after* their deterministic
    //descendants (i.e. in reverse topological order).
    list<DeterministicNode*> const *dc = node->deterministicChildren();
    for(list<DeterministicNode*>::const_iterator p = dc->begin();
	p != dc->end(); ++p)
    {
	if (hasDescendant(*p, target_set, dnodes, known_dnodes)) {
	    ans = true;
	}
    }

    //At this point we have visited all deterministic descendants and can
    //classify the current node.
    known_dnodes[node] = ans;
    if (ans) {
	dnodes.push_back(node);
    }
    return ans;
}


/* hasDescendant for stochastic nodes */
static bool hasDescendant(StochasticNode *node, 
			  set<Node const*> const &target_set,
			  vector<DeterministicNode*> &dnodes,
			  map<DeterministicNode*, bool> &known_dnodes)
{
    //This is a stripped-down version of hasDescendant for deterministic
    //nodes. We do not keep track of known stochastic nodes
    bool ans = target_set.count(node);

    list<DeterministicNode*> const *dc = node->deterministicChildren();
    for(list<DeterministicNode*>::const_iterator p = dc->begin();
	p != dc->end(); ++p)
    {
	if (hasDescendant(*p, target_set, dnodes, known_dnodes)) {
	    ans = true;
	}
    }

    return ans;
}

/*  Find stochastic parents of given set of nodes within the model,
    and their intermediate deterministic descendants in topological
    order */
static void findStochasticParents(vector<Node const *> const &tgt_nodes, 
				  Model const &model,
				  vector<StochasticNode *> &stoch_parents,
				  vector<DeterministicNode *> &dtrm_parents)
{
    //Copy input vector tgt_nodes to a set
    set<Node const *> tgt_set;
    for (unsigned int i = 0; i < tgt_nodes.size(); ++i) {
	tgt_set.insert(tgt_nodes[i]);
    }

    //Set up auxliary variables required for hasDescendant and call
    map<DeterministicNode*, bool> known_dnodes;
    vector<StochasticNode*> const &snodes = model.stochasticNodes();
    for (unsigned int i = 0; i < snodes.size(); ++i) {
	if (hasDescendant(snodes[i], tgt_set, dtrm_parents, known_dnodes)) {
	    stoch_parents.push_back(snodes[i]);
	}
    }
    
    // Nodes are pushed onto dtrm_parents in reverse topological order.
    // Reverse the vector to get them in topological (forward sampling) order
    reverse(dtrm_parents.begin(), dtrm_parents.end());
}

    /*
static void cloneNodes(vector<StochasticNode*> const &nodes, 
		       vector<StochasticNode*> &newnodes,
		       map<Node const*, Node const*> &remap)
{
    // Clone nodes without remapping parents
    for (unsigned int i = 0; i < nodes.size(); ++i) {
	StochasticNode *newnode = nodes[i]->clone(nodes[i]->parents());
	newnodes.push_back(newnode);
	remap[nodes[i]] = newnode;
    }
}

static void cloneNodes(vector<DeterministicNode*> const &nodes, 
		       vector<DeterministicNode*> &newnodes,
		       map<Node const*, Node const*> &remap)
{
    // Clone nodes with remapped parents
    for (unsigned int i = 0; i < nodes.size(); ++i) {
	vector<Node const *> args = nodes[i]->parents();
	for (unsigned int j = 0; j < args.size(); ++j) {
	    map<Node const*, Node const*>::const_iterator p =
		remap.find(args[j]);
	    if (p != remap.end()) 
		args[j] = p->second;
	}
	DeterministicNode *newnode = nodes[i]->clone(args);
	newnodes.push_back(newnode);
	remap[nodes[i]] = newnode;
    }
}
    */

static Node* 
getMixtureNode1(NodeArray *array, vector<SSI> const &limits, Compiler *compiler)
// Try to simplify the mixture node by enumerating all possible values
// that the indices can take. This can only be done under certain
// conditions.
{
    unsigned int ndim = limits.size();

    vector<Node const*> indices;
    for (unsigned int i = 0; i < ndim; ++i) {
	if (limits[i].node) {
	    indices.push_back(limits[i].node);
	}
    }
    unsigned int nvi = indices.size();

    vector<StochasticNode *> sparents;
    vector<DeterministicNode *> dparents;
    findStochasticParents(indices, compiler->model(), sparents, dparents);

    unsigned int nparents = sparents.size();
    if (nparents > 10) {
	//This algorithm grinds to a halt with too many stochastic
	//parents. So bail out after 10
	return 0;
    }
    vector<int> lower(nparents), upper(nparents);
    for (unsigned int i = 0; i < nparents; ++i) {
	StochasticNode const *snode = sparents[i];

	if (snode->length() != 1 || !snode->isDiscreteValued() ||
	    !isSupportFixed(snode))
	{
	    //Must have discrete, scalar parents with fixed support
	    return 0;
	}
	
	// Get lower and upper limits of support
	double l = JAGS_NEGINF, u = JAGS_POSINF;
	snode->support(&l, &u, 1U, 0);
	if (!jags_finite(l) || !jags_finite(u)) {
	    return 0; //Unbounded parent => serious trouble
	}
	if (l < -INT_MAX || u > INT_MAX) {
	    return 0; //Can't cast to int
	}
	lower[i] = static_cast<int>(l);
	upper[i] = static_cast<int>(u);
    }

    //Store current value of all stochastic parents and deterministic nodes
    vector<vector<double> > stored_parent_values, stored_dtrm_values;
    save(sparents, stored_parent_values);
    save(dparents, stored_dtrm_values);

    // Create a set containing all possible values that the stochastic
    // indices can take
    set<vector<int> > index_values;
    vector<int>  this_index(indices.size(),1);
    SimpleRange stoch_node_range(lower, upper);

    for (RangeIterator i(stoch_node_range); !i.atEnd(); i.nextLeft()) {

	for (unsigned int j = 0; j < sparents.size(); ++j) {
	    double v = i[j];
	    sparents[j]->setValue(&v, 1, 0);
	}

	for (unsigned int k = 0; k < dparents.size(); ++k) {
	    dparents[k]->deterministicSample(0);
	}
    
	for (unsigned int l = 0; l < indices.size(); ++l) {
	    this_index[l] = static_cast<int>(*indices[l]->value(0));
	}
    
	index_values.insert(this_index);
    }

    //Restore value of all stochastic parents and deterministic nodes
    restore(sparents, stored_parent_values);
    restore(dparents, stored_dtrm_values);
    
    // Now set up the possible subsets defined by the stochastic indices
    vector<int> variable_offset;
    vector<vector<int> > scope(ndim);
    for (unsigned int j = 0; j < ndim; ++j) {
	if (limits[j].node) {
	    variable_offset.push_back(j);
	}
	else {
	    scope[j] = limits[j].indices;
	}
    }

    vector<pair<vector<int>, Range> > ranges;  
    set<vector<int> >::const_iterator p;
    for (p = index_values.begin(); p != index_values.end(); ++p) {

	vector<int> const &i = *p;
	for (unsigned int k = 0; k < nvi; ++k) {
	    scope[variable_offset[k]] = vector<int>(1, i[k]);
	}
	ranges.push_back(pair<vector<int>, Range>(i, Range(scope)));
    }

    //Look out for trivial mixture nodes in which all subsets are the same.
    bool trivial = true;
    Node *subset_node0 = array->getSubset(ranges[0].second, compiler->model());

    map<vector<int>, Node const *> subsets;  
    for (unsigned int i = 0; i < ranges.size(); ++i) {
	Node *subset_node = array->getSubset(ranges[i].second, 
					     compiler->model());
	if (!subset_node) {
	    return 0;
	}
	subsets[ranges[i].first] = subset_node;
	if (subset_node != subset_node0) {
	    trivial = false;
	}
    }

    if (trivial) {
	// Nothing to do here! All subset nodes are the same, so we
	// just return the subset node.
	return subset_node0;
    }

    return compiler->mixtureFactory1().getMixtureNode(indices, subsets, 
						      compiler->model());
}

static Node * 
getMixtureNode2(NodeArray *array, vector<SSI> const &limits, Compiler *compiler)
{
    vector<pair<vector<int>, Range> > ranges;  
    getSubsetRanges(ranges, limits, array->range());

    map<vector<int>, Node const *> mixmap;
    for (unsigned int i = 0; i < ranges.size(); ++i) {
	Node *subset_node =
	    array->getSubset(ranges[i].second, compiler->model());
	if (subset_node) {
	    mixmap[ranges[i].first] = subset_node;
	}
	else {
	    return 0;
	}

    }

    vector<Node const *> indices;
    for (unsigned int i = 0; i < limits.size(); ++i) {
	if (limits[i].node) {
	    indices.push_back(limits[i].node);
	}
    }
    
    return compiler->mixtureFactory2().getMixtureNode(indices, mixmap, 
						      compiler->model());
}

Node * getMixtureNode(ParseTree const * var, Compiler *compiler)
{
    if (var->treeClass() != P_VAR) {
	throw logic_error("Expecting variable expression");
    }
  
    NodeArray *array = compiler->model().symtab().getVariable(var->name());
    if (array == 0) {
	throw runtime_error(string("Unknown parameter: ") + var->name());
    }

    vector<ParseTree*> const &range_list = var->parameters();
    vector<SSI> limits;
    unsigned int ndim = array->range().ndim(false);
    if (range_list.size() != ndim) {
	throw runtime_error("Dimension mismatch taking variable subset of " + 
			    var->name());
    }

    unsigned int nvi = 0; //Count number of variable indices
    for (unsigned int i = 0; i < ndim; ++i) {
	ParseTree const *range_element = range_list[i];
	if (range_element->treeClass() != P_RANGE) {
	    throw runtime_error("Malformed range expression");
	}
	SSI ssi;
	ParseTree const *p0;
	switch(range_element->parameters().size()) {
	case 0:
	    // Index is empty, implying the whole range 
	    ssi.indices = array->range().scope()[i];
	    break;
	case 1:
	    p0 = range_element->parameters()[0];
	    if(!compiler->indexExpression(p0, ssi.indices)) {
		ssi.node = compiler->getParameter(p0);
		if (ssi.node == 0)
		    return 0;
		else
		    ++nvi;
	    }
	    break;
	default:
	    throw logic_error("Invalid range expression");
	}
	//Check validity of subset index
	if (ssi.node) {
	    if (!ssi.node->isDiscreteValued()) {
		throw NodeError(ssi.node, "Continuous node used as index");
	    }
	    if (ssi.node->length() != 1) {
		throw NodeError(ssi.node, "Vector node used as index");
	    }
	}
	else {
	    if (ssi.indices.empty()) {
		throw runtime_error("Requested invalid variable subset of " +
				    var->name());
	    }
	    for (unsigned int j = 0; j < ssi.indices.size(); ++j) {
		if (ssi.indices[j] < array->range().lower()[i] ||
		    ssi.indices[j] > array->range().upper()[i])
		{
		    throw runtime_error("Requested invalid variable subset of "
					+ var->name());
		}
	    }
	}
	limits.push_back(ssi);
    }
    
    //Check number of variable indices (nvi)
    if (nvi == 0) {
	throw logic_error("Trivial mixture node");
    }

    Node *mnode = getMixtureNode1(array, limits, compiler);
    if (mnode)
	return mnode;
    else
	return getMixtureNode2(array, limits, compiler);
}


    /**
     * This is a diagnostic function that is called by the compiler
     * in _resolution_level 2 when a mixture node cannot be resolved.
     *
     * FIXME: There is some serious cut-and-paste in this function.
     * However, this is currently the simplest way to implement these
     * checks without perturbing the API
     */
    void getMissingMixParams(ParseTree const * var, UMap &umap,
			     Compiler *compiler)
    {
	NodeArray *array = compiler->model().symtab().getVariable(var->name());

	vector<ParseTree*> const &range_list = var->parameters();
	vector<SSI> limits;
	unsigned int ndim = array->range().ndim(false);

	bool resolved = true;
	for (unsigned int i = 0; i < ndim; ++i) {
	    ParseTree const *range_element = range_list[i];
	    SSI ssi;
	    ParseTree const *p0;
	    switch(range_element->parameters().size()) {
	    case 0:
		// Index is empty, implying the whole range 
		ssi.indices = array->range().scope()[i];
		break;
	    case 1:
		p0 = range_element->parameters()[0];
		if(!compiler->indexExpression(p0, ssi.indices)) {
		    ssi.node = compiler->getParameter(p0);
		    if (ssi.node == 0) {
			resolved = false;
		    }
		}
		break;
	    default:
		throw logic_error("Invalid range expression");
	    }

	    limits.push_back(ssi);
	}

	if (!resolved) return;

	vector<pair<vector<int>, Range> > ranges;  
	getSubsetRanges(ranges, limits, array->range());

	for (unsigned int i = 0; i < ranges.size(); ++i) {
	    Range const &subset_range = ranges[i].second;
	    Node *node = array->getSubset(subset_range, compiler->model());
	    if (node) continue;

	    /* Make a note of all subsets that could not be resolved. */
	    bool empty = true;
	    for (RangeIterator r(subset_range); !r.atEnd(); r.nextLeft()) {
		if (array->getSubset(SimpleRange(r,r), compiler->model())) {
		    empty = false;
		    break;
		}
	    }
	    if (empty) {
		//No elements found, so report the whole range
		//as missing, e.g. x[1:5]
		pair<string, Range> upair(var->name(), subset_range);
		if (umap.find(upair) == umap.end()) {
		    set<int> lines;
		    lines.insert(var->line());
		    umap[upair] = lines;
		}
		else {
		    umap.find(upair)->second.insert(var->line());
		}
	    }
	    else {
		//Some elements found, so report only the missing
		//elements, e.g. x[5]
		for (RangeIterator r(subset_range); !r.atEnd(); r.nextLeft()) {
		    SimpleRange sr(r,r);
		    if (!array->getSubset(sr, compiler->model())) {
			pair<string, Range> upair(var->name(), sr);
			if (umap.find(upair) == umap.end()) {
			    set<int> lines;
			    lines.insert(var->line());
			    umap[upair] = lines;
			}
			else {
			    umap.find(upair)->second.insert(var->line());
			}
		    }
		}
	    }
	}
    }    


} //namespace jags
