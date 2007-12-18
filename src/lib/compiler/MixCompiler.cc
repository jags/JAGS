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
   the case where all the direct Stochastic ancestors of Y are discrete
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
#include <sampler/Sampler.h>
#include <util/nainf.h>

#include "MixCompiler.h"

#include <stdexcept>
#include <string>
#include <set>
#include <climits>

using std::vector;
using std::pair;
using std::logic_error;
using std::runtime_error;
using std::string;
using std::set;
using std::map;

//Structure to hold subset indices
struct SSI {
  Node const *node;
  int lower;
  int upper;
  SSI() : node(0), lower(0), upper(0) {};
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
			    Range const &default_range)
{
  unsigned int ndim = limits.size();
  //Count number of variable indices (nvi)
  int nvi = 0;
  for (unsigned int j = 0; j < ndim; ++j) {
    if (limits[j].node != 0)
      ++nvi;
  }

  // Create upper and lower bounds
  vector<int> variable_offset(nvi,1), variable_lower(nvi,1), variable_upper(nvi,1);
  vector<int> lower_index(ndim,1), upper_index(ndim,1);
  int k = 0;
  for (unsigned int j = 0; j < ndim; ++j) {
    if (limits[j].node != 0) {
      variable_offset[k] = j;
      variable_lower[k] = default_range.lower()[j];
      variable_upper[k] = default_range.upper()[j];
      ++k;
    }
    else {
      lower_index[j] = limits[j].lower;
      upper_index[j] = limits[j].upper;
    }
  }
  
  for (RangeIterator i(Range(variable_lower, variable_upper)); !i.atEnd(); 
       i.nextLeft()) {
    for (int k = 0; k < nvi; ++k) {
      lower_index[variable_offset[k]] = i[k];
      upper_index[variable_offset[k]] = i[k];
    }
    subsets.push_back(pair<vector<int>, Range>(i, Range(lower_index, upper_index)));
  }
}

static void markParents(Node const *node, GraphMarks &marks)
{
    /* Recursively mark node and its parents until we get back to the
       nearest stochastic ancestors */

    if (node->isRandomVariable()) {
        if (asStochastic(node)) 
	    marks.mark(node, 2);
    }
    else if (marks.mark(node) == 0) {
	marks.mark(node, 1);
	vector<Node const*>::const_iterator p;
	for (p = node->parents().begin(); p != node->parents().end(); ++p) {
	    markParents(*p, marks);
	}
    }
}

static void classifyParents(vector<Node const *> const &indices, 
			    Graph const &graph,
			    vector<StochasticNode *> &stoch_nodes,
			    vector<Node *> &dtrm_nodes)
{
    GraphMarks marks(graph);

    /* Classify parents of each node */
    vector<Node const*>::const_iterator p;
    for (p = indices.begin(); p != indices.end(); ++p) {
	markParents(*p, marks);
    }
    
    Graph dgraph;
    set<Node *>::const_iterator q;
    for (q = graph.nodes().begin(); q != graph.nodes().end(); ++q) {
	switch(marks.mark(*q)) {
	case 0:
	    break;
	case 1:
	    dgraph.add(*q);
	    break;
	case 2:
	    stoch_nodes.push_back(static_cast<StochasticNode*>(*q));	    
	    break;
	default:
	    break;
	}
    }
    /* Get deterministic nodes in forward sampling order */
    dgraph.getSortedNodes(dtrm_nodes);
}

static Node* 
getMixtureNode1(NodeArray *array, vector<SSI> const &limits, Compiler *compiler)
/* Try to simplify the mixture node by enumerating all possible values
   that the indices can take. This can only be done under certain
   conditions. */
{
    unsigned int ndim = limits.size();

    vector<Node const*> indices;
    unsigned int nvi = 0;
    for (unsigned int i = 0; i < ndim; ++i) {
	if (limits[i].node) {
	    indices.push_back(limits[i].node);
	    ++nvi;
	}
    }

    Graph &compiler_graph = compiler->model().graph();

    /* Find stochastic parents (ancestors) and deterministic ancestors in
       forward sampling order */
  
    vector<StochasticNode*> stoch_parents;
    vector<Node*> dtrm_nodes;
    classifyParents(indices, compiler_graph, stoch_parents, dtrm_nodes);

    /* Test to see if all stochastic parents are discrete, scalar, with
       bounded support. If not we give up. */

    unsigned int nparents = stoch_parents.size();  
    for (unsigned int i = 0; i < nparents; ++i) {
	StochasticNode const *snode = stoch_parents[i];
	if (snode->length() != 1 || !snode->isDiscreteValued() ||
	    snode->isObserved() || isBounded(snode)) 
	{
	    return 0;
	}
    }

    vector<int> lower(nparents,1), upper(nparents,1);
    for (unsigned int i = 0; i < nparents; ++i) {
	StochasticNode const *snode = stoch_parents[i];

	/* Check that support of node is fixed */
	if (!isSupportFixed(snode)) {
	    return 0;
	}

	/* To be safe, we cycle over all chains */
	for (unsigned int n = 0; n < snode->nchain(); ++n) {
	    // Get lower and upper limits of support
	    double l = JAGS_NEGINF, u = JAGS_POSINF;
	    support(&l, &u, 1U, snode, n);
	    if (!jags_finite(l) || !jags_finite(u)) {
		return 0; //Unbounded parent => serious trouble
	    }
	    if (l < -INT_MAX || u > INT_MAX) {
		return 0; //Can't cast to int
	    }
	    int il = static_cast<int>(l);
	    int iu = static_cast<int>(u);
	    if (n == 0 || il < lower[i]) {
		lower[i] = il;
	    }
	    if (n == 0 || iu > upper[i]) {
		upper[i] = iu;
	    }
	}
    }

    /* Save current stochastic node values in chain 0 */
    vector<double> vsave(stoch_parents.size());
    for (unsigned int j = 0; j < stoch_parents.size(); ++j) {
	vsave[j] = *stoch_parents[j]->value(0);
    }

    /* Create a set containing all possible values that the stochastic
       indices can take */

    set<vector<int> > index_values;
    vector<int>  this_index(indices.size(),1);
    Range stoch_node_range(lower, upper);

    for (RangeIterator i(stoch_node_range); !i.atEnd(); i.nextLeft()) {

	for (unsigned int j = 0; j < stoch_parents.size(); ++j) {
	    double v = i[j];
	    stoch_parents[j]->setValue(&v, 1, 0);
	}

	for (unsigned int k = 0; k < dtrm_nodes.size(); ++k) {
	    dtrm_nodes[k]->deterministicSample(0);
	}
    
	for (unsigned int l = 0; l < indices.size(); ++l) {
	    this_index[l] = static_cast<int>(*indices[l]->value(0));
	}
    
	index_values.insert(this_index);
    }

    /* Restore saved values to chain 0 */
    for (unsigned int j = 0; j < stoch_parents.size(); ++j) {
	double v = vsave[j];
	stoch_parents[j]->setValue(&v, 1, 0);
    }
    for (unsigned int k = 0; k < dtrm_nodes.size(); ++k) {
	dtrm_nodes[k]->deterministicSample(0);
    }

    /* Now set up the possible subsets defined by the stochastic indices */

    vector<int> variable_offset(nvi, 1);
    vector<int> lower_index(ndim, 1), upper_index(ndim, 1);
    int k = 0;
    for (unsigned int j = 0; j < ndim; ++j) {
	if (limits[j].node != 0) {
	    variable_offset[k++] = j;
	}
	else {
	    lower_index[j] = limits[j].lower;
	    upper_index[j] = limits[j].upper;
	}
    }

    vector<pair<vector<int>, Range> > ranges;  
    set<vector<int> >::const_iterator p;
    for (p = index_values.begin(); p != index_values.end(); ++p) {

	vector<int> const &i = *p;
	for (unsigned int k = 0; k < nvi; ++k) {
	    lower_index[variable_offset[k]] = i[k];
	    upper_index[variable_offset[k]] = i[k];
	}
	ranges.push_back(pair<vector<int>, Range>(i, Range(lower_index, upper_index)));
    }

    /* Convert these into subsets */

    //Look out for trivial mixture nodes in which all subsets are the same.
    bool trivial = true;
    Node *subset_node0 = array->getSubset(ranges[0].second, compiler_graph);

    map<vector<int>, Node const *> subsets;  
    for (unsigned int i = 0; i < ranges.size(); ++i) {
	Node *subset_node = array->getSubset(ranges[i].second, compiler_graph);
	if (!subset_node)
	    return 0;
	subsets[ranges[i].first] = subset_node;
	if (subset_node != subset_node0)
	    trivial = false;
    }

    if (trivial) {
	// Nothing to do here! All subset nodes are the same, so we
	// just return the subset node.
	return subset_node0;
    }

    return compiler->mixtureFactory().getMixtureNode(indices, subsets, 
						     compiler_graph);
}

static Node * 
getMixtureNode2(NodeArray *array, vector<SSI> const &limits, Compiler *compiler)
{
    Graph &cgraph = compiler->model().graph();

    vector<pair<vector<int>, Range> > ranges;  
    getSubsetRanges(ranges, limits, array->range());

    map<vector<int>, Node const *> subsets;
    for (unsigned int i = 0; i < ranges.size(); ++i) {
	Node *subset_node =   array->getSubset(ranges[i].second, cgraph);
	if (!subset_node)
	    return 0;
	subsets[ranges[i].first] = subset_node;
    }

    vector<Node const *> indices;
    for (unsigned int i = 0; i < limits.size(); ++i) {
	if (limits[i].node) {
	    indices.push_back(limits[i].node);
	}
    }

    return compiler->mixtureFactory().getMixtureNode(indices, subsets, cgraph);
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
    ssi.node = 0;
    ParseTree const *p0, *p1;
    switch(range_element->parameters().size()) {
    case 0:
      // Index is empty, implying the whole range 
      ssi.lower = array->range().lower()[i];
      ssi.upper = array->range().upper()[i];
      break;
    case 1:
      // Single index: upper = lower
      p0 = range_element->parameters()[0];
      if(compiler->indexExpression(p0, ssi.lower)) {
	ssi.upper = ssi.lower;
      }
      else {
	  ssi.node = compiler->getParameter(p0);
	if (ssi.node == 0)
	  return 0;
	else
	  ++nvi;
      }
      break;
    case 2:
      // Upper and lower indices
      p0 = range_element->parameters()[0];
      p1 = range_element->parameters()[1];
      if(compiler->indexExpression(p0, ssi.lower)) {
	if (!compiler->indexExpression(p1, ssi.upper)) {
	  return 0;
	}
      }
      else {
	  ssi.node = compiler->getParameter(p0);
	  if (compiler->getParameter(p1) != ssi.node)
	  return 0;
	else
	  ++nvi;
      }
      break;
    default:
      throw logic_error("Invalid range expression");
    }
    //Check validity of limits
    if (ssi.node == 0) {
      if (ssi.lower < array->range().lower()[i] ||
	  ssi.upper > array->range().upper()[i] ||
	  ssi.upper < ssi.lower)
	{
	  throw runtime_error("Requested invalid variable subset of " + 
			      var->name());
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
