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
#include <graph/DeterministicNode.h>
#include <distribution/DistScalar.h>
#include <model/NodeArray.h>
#include <graph/MixtureNode.h>

#include "MixCompiler.h"

#include <stdexcept>
#include <string>
#include <set>

using std::vector;
using std::pair;
using std::logic_error;
using std::runtime_error;
using std::string;
using std::set;

//Structure to hold subset indices
struct SSI {
  Node *node;
  long lower;
  long upper;
};


/* 
   Given a vector of subset indices, this function modifies the argument
   "subsets", so that it contains a vector of pairs:
   - Index: Possible value of the variable indices
   - Range: Corresponding range of the array to take as a subset

   All possible values of Index are included in the vector, based on the
   the range of the array we are taking subsets from (default_range).
*/
static void getSubsetRanges(vector<pair<Index, Range> > &subsets,  
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
  Index variable_offset(nvi), variable_lower(nvi), variable_upper(nvi);
  Index lower_index(ndim), upper_index(ndim);
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
    subsets.push_back(pair<Index, Range>(i, Range(lower_index, upper_index)));
  }
}

/* Helper function for classifyParents */
static void classifyNode(Node *node, Graph &sgraph, Graph &dgraph)
{
    if (node->isVariable()) {
        if (asStochastic(node)) 
	   sgraph.add(node);
    }
    else if (!dgraph.contains(node)) {
	dgraph.add(node);      
	for (vector<Node*>::const_iterator p = node->parents().begin();
	     p != node->parents().end(); ++p) {
	    classifyNode(*p, sgraph, dgraph);
	}
    }
}

/* 
   Trace ancestors of Stochastic nodes in the vector "indices" back to
   their stochastic ancestors. This function works exactly like
   Sampler::classifyChildren, but for parents rather than children.
*/
static void classifyParents(vector<Node *> const &indices,
			    vector<StochasticNode *> &stoch_nodes,
			    vector<DeterministicNode *> &dtrm_nodes)
{
  Graph dgraph, sgraph;

  /* Classify parents of each node */
  for (vector<Node*>::const_iterator p = indices.begin(); p != indices.end(); ++p) {
    classifyNode(*p, sgraph, dgraph);
  }
  
  vector<Node*> svector;
  sgraph.getNodes(svector);
  for (vector<Node*>::iterator i = svector.begin(); i != svector.end(); 
       ++i) 
    {
      //FIXME: downcasting non-constant pointer
      stoch_nodes.push_back(static_cast<StochasticNode*>(*i));
    }
  
  vector<Node*> dvector;
  dgraph.getSortedNodes(dvector);
  for (vector<Node*>::iterator i = dvector.begin(); i != dvector.end();
       ++i)
    {
      dtrm_nodes.push_back(dynamic_cast<DeterministicNode*>(*i));
    }
}

static Node* getMixtureNode1(NodeArray *array, vector<SSI> const &limits, Compiler *compiler)
/* Try to simplify the mixture node by enumerating all possible values
   that the indices can take. This can only be done under certain
   conditions. */
{
  unsigned int ndim = limits.size();

  vector<Node*> indices;
  unsigned int nvi = 0;
  for (unsigned int i = 0; i < ndim; ++i) {
    if (limits[i].node) {
      indices.push_back(limits[i].node);
      ++nvi;
    }
  }
  
  /* Find stochastic parents (ancestors) and deterministic ancestors in
     forward sampling order */
  
  vector<StochasticNode*> stoch_parents;
  vector<DeterministicNode*> dtrm_nodes;
  classifyParents(indices, stoch_parents, dtrm_nodes);
  
  /* Test to see if all stochastic parents are discrete, scalar, with
     bounded support. If not we give up. */
  
  bool cansimplify = true;
  for (unsigned int i = 0; i < stoch_parents.size(); ++i) {
    StochasticNode const *snode = stoch_parents[i];
    if (snode->length() != 1 || !snode->isDiscreteValued() ||
	snode->isObserved() || snode->isBounded()) 
      {
	cansimplify = false;
	break;
      }
    if (dynamic_cast<DistScalar const *>(snode->distribution()) == 0)
      {
	/* We need to be able to downcast the distribution to DistScalar
	   in order to access the "l" and "u" member functions */
	cansimplify = false;
	break;
      }
  }
  if (!cansimplify) {
    return 0; 
  }

  unsigned long nparents = stoch_parents.size();
  Index lower(nparents), upper(nparents);
  for (unsigned int i = 0; i < nparents; ++i) {
    StochasticNode const *snode = stoch_parents[i];
    Distribution const *dist = snode->distribution();
    DistScalar const *dscalar = dynamic_cast<DistScalar const *>(dist);

    /* To be safe, we cycle over all chains */
    for (unsigned int n = 0; n < snode->nchain(); ++n) {
      // Get absolute lower and upper bounds
      double l = dscalar->l(snode->parameters(n), true);
      double u = dscalar->u(snode->parameters(n), true);
      if (l == -DBL_MAX || u == DBL_MAX) {
	return 0; //Unbounded parent => serious trouble
      }
      if (l < -LONG_MAX || u > LONG_MAX) {
	return 0; //Can't cast to long
      }
      long il = static_cast<long>(l);
      long iu = static_cast<long>(u);
      if (n == 0 || il < lower[i]) {
	lower[i] = il;
      }
      if (n == 0 || iu > upper[i]) {
	upper[i] = iu;
      }
    }
  }

  //Save current node value in chain 0;
  vector<double> vsave(stoch_parents.size());
  for (unsigned int j = 0; j < stoch_parents.size(); ++j) {
    vsave[j] = *stoch_parents[j]->value(0);
  }

  /* Create a set containing all possible values that the stochastic
     indices can take */

  set<Index> index_values;
  Index this_index(indices.size());
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
      this_index[l] = static_cast<long>(*indices[l]->value(0));
    }
    
    index_values.insert(this_index);
  }

  /* Restore saved values to chain 0 */
  for (unsigned int j = 0; j < stoch_parents.size(); ++j) {
    double v = vsave[j];
    stoch_parents[j]->setValue(&v, 1, 0);
  }

  /* Now set up the possible subsets defined by the stochastic indices */

  Index variable_offset(nvi);
  Index lower_index(ndim), upper_index(ndim);
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

  vector<pair<Index, Range> > ranges;  
  set<Index>::const_iterator p;
  for (p = index_values.begin(); p != index_values.end(); ++p) {

    Index const &i = *p;
    for (unsigned int k = 0; k < nvi; ++k) {
      lower_index[variable_offset[k]] = i[k];
      upper_index[variable_offset[k]] = i[k];
    }
    ranges.push_back(pair<Index, Range>(i, Range(lower_index, upper_index)));
  }

  /* Convert these into subsets */

  vector<pair<Index, Node*> > subsets;  
  for (unsigned int i = 0; i < ranges.size(); ++i) {
    Node *subset_node = array->getSubset(ranges[i].second);
    if (!subset_node)
      return 0;
    subsets.push_back(pair<Index, Node*>(ranges[i].first, subset_node));
  }

  return compiler->mixtureFactory().getMixtureNode(indices, subsets);
}

static Node * getMixtureNode2(NodeArray *array, vector<SSI> const &limits, Compiler *compiler)
{
  vector<pair<Index, Range> > ranges;  
  getSubsetRanges(ranges, limits, array->range());
  vector<pair<Index, Node*> > subsets;  
  for (unsigned int i = 0; i < ranges.size(); ++i) {
    Node *subset_node = array->getSubset(ranges[i].second);
    if (!subset_node)
      return 0;
    subsets.push_back(pair<Index, Node*>(ranges[i].first, subset_node));
  }

  vector<Node*> indices;
  for (unsigned int i = 0; i < limits.size(); ++i) {
    if (limits[i].node) {
      indices.push_back(limits[i].node);
    }
  }

  return compiler->mixtureFactory().getMixtureNode(indices, subsets);
}

Node * getMixtureNode(ParseTree const * var, Compiler *compiler)
{
  if (var->treeClass() != P_VAR) {
    throw logic_error("Expecting variable expression");
  }
  
  NodeArray *array = compiler->symTab().getVariable(var->name());
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
