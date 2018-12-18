#include <config.h>
#include <graph/AggNode.h>
#include <graph/GraphMarks.h>
#include <graph/Graph.h>

#include <vector>
#include <stdexcept>

using std::vector;
using std::set;
using std::logic_error;
using std::length_error;
using std::out_of_range;
using std::string;

namespace jags {

/* 
   Two utility functions "substitute_parents" and "substitute_offsets"
   ensure that an AggNode is never the parent of another AggNode, even
   if it was originally defined this way. The substitution makes
   subsetting of aggregate nodes more efficient. For example:
   
   X <- A[4:8]
   Y <- X[3:4]
   
   is equivalent to
   
   X <- A[4:8]
   Y <- A[7:8]
   
   Without the substitute_* functions, node A would be a grand-parent
   of node Y instead of a parent.  
*/

static vector<Node const *> sub_parents(vector<Node const *> const &parents, 
					vector<unsigned long> const &offsets)
{
    // Substitute parent nodes that are themselves AggNodes
    vector<Node const *> newparents(parents);
    for (unsigned long i = 0; i < parents.size(); i++) {
	AggNode const *aggpar = dynamic_cast<AggNode const *>(parents[i]);
	if (aggpar) {
	    newparents[i] = aggpar->parents()[offsets[i]];
	}
    }
    return newparents;
}

static vector<unsigned long> sub_offsets(vector<Node const *> const &parents, 
					vector<unsigned long> const &offsets)
{
    // Substitute offsets for parents that are themselves AggNodes

    vector<unsigned long> newoffsets(offsets);
    for (unsigned long i = 0; i < offsets.size(); i++) {
	AggNode const *aggpar = dynamic_cast<AggNode const *>(parents[i]);
	if (aggpar) {
	    newoffsets[i] = aggpar->offsets()[offsets[i]];
	}
    }
    return newoffsets;
}

AggNode::AggNode(vector<unsigned long> const &dim,
		 unsigned int nchain,
		 vector<Node const *> const &parents,
                 vector<unsigned long> const &offsets)
    : DeterministicNode(dim, nchain, sub_parents(parents, offsets)), 
      _offsets(sub_offsets(parents, offsets)),
      _parent_values(_length * nchain), _discrete(true)
{
    // Check argument lengths
    if (_length != parents.size() || _length != offsets.size()) {
	throw length_error ("Length mismatch in Aggregate Node constructor");
    }

    /* 
       Note that we cannot use the original arguments "parents" and "offsets"
       due to possible substitution. Use this->parents() and _offsets instead.
    */
    vector<Node const *> const &par = this->parents();

    // Check that offsets are valid
    for (unsigned long i = 0; i < _length; i++) {
	if (_offsets[i] >= par[i]->length())
	    throw out_of_range("Invalid offset in Aggregate Node constructor");
    }
  
    // Setup parent values
    for (unsigned int ch = 0; ch < _nchain; ++ch) {
	for (unsigned long i = 0; i < _length; ++i) {
	    _parent_values[i + ch * _length] = par[i]->value(ch) + _offsets[i];
	}
    }

    // Check discreteness
    for (unsigned long i = 0; i < par.size(); ++i) {
	if (!par[i]->isDiscreteValued()) {
	    _discrete = false;
	    break;
	}
    }

    // Initialize if fixed. See comments in the DeterministicNode constructor
    if (isFixed()) {
	for (unsigned int ch = 0; ch < _nchain; ++ch) {
	    deterministicSample(ch);
	}
    }

}

AggNode::~AggNode()
{
}

void AggNode::deterministicSample(unsigned int chain)
{
    unsigned long N = _length * chain;
    for (unsigned long i = 0; i < _length; ++i) {
	_data[i + N] = *_parent_values[i + N];
    }
}

    bool AggNode::isDifferentiable(Node const *arg) const
    {
	return true;
    }
    
    void AggNode::gradient(double *grad, Node const *arg,
			   unsigned int chain) const
    {
	//FIXME put in some check code here...
	
	unsigned long P = _length;

	auto par = parents();
	for (unsigned int p = 0; p < P; ++p) {
	    if (par[p] == arg) {
		unsigned long q = _offsets[p];
		grad[p + P * q] += 1;
	    }
	}
    }

    
/*
bool AggNode::isLinear(GraphMarks const &linear_marks, bool fixed) const
{
    vector<Node const *> const &par = parents();
    for (unsigned int i = 0; i < par.size(); ++i) {
	if (linear_marks.graph().contains(par[i])) {
	    switch(linear_marks.mark(par[i])) {
	    case MARK_NULL: case MARK_TRUE:
		break;
	    case MARK_FALSE:
		return false;
		break;
	    default:
		throw logic_error("Invalid graphmarks in AggNode::isLinear");
	    }
	}
    }
    return true;
}
*/

bool AggNode::isClosed(set<Node const *> const &ancestors, 
		       ClosedFuncClass fc, bool fixed) const
{
    
    switch(fc) {
    case DNODE_ADDITIVE:
    {
	//Only one parent may be additive, and it must be completely
	//embedded in the AggNode
	Node const *pnode = nullptr;
	vector<bool> pmask;
	vector<Node const *> const &par = parents();
	for (unsigned int i = 0; i < par.size(); ++i) {
	    if (ancestors.count(par[i])) {
		if (pnode == nullptr) {
		    pnode = par[i];
		    pmask = vector<bool>(pnode->length(), false);
		}
		else {
		    if (par[i] != pnode) return false;
		    if (pmask[_offsets[i]]) return false;
		}
		pmask[_offsets[i]] = true;
	    }
	    else if (fixed) {
		if (!par[i]->isFixed()) return false;
	    }
	}
	if (!allTrue(pmask)) return false;
    }
    break;
    case DNODE_SCALE:
	//All parents must be scale transformations
	for (unsigned int i = 0; i < parents().size(); ++i) {
	    if (ancestors.count(parents()[i]) == 0) {
		return false;
	    }
	}
	break;
    case DNODE_SCALE_MIX:
	//The aggnode must be a subset so only one distinct parent
	if (ancestors.count(parents()[0]) == 0)
	    return false;
	for (unsigned int i = 1; i < parents().size(); ++i) {
	    if (parents()[i] != parents()[0])
		return false;
	}
	break;
    case DNODE_LINEAR:
	break;
    case DNODE_POWER:
	return false;
    }
    return true;
}

bool AggNode::checkParentValues(unsigned int) const
{
   return true;
}

string AggNode::deparse(vector<string> const &parents) const
{
    return string("aggregate(") + parents.front() + "..." + 
	parents.back() + ")";
}

/*    
DeterministicNode *
AggNode::clone(vector<Node const *> const &parents) const
{
    return new AggNode(this->dim(), parents, _offsets);
}
*/

bool AggNode::isDiscreteValued() const
{
  return _discrete;
}

vector<unsigned long> const &AggNode::offsets() const
{
    return _offsets;
}

} //namespace jags
