#include <config.h>
#include <graph/DeterministicNode.h>

using std::vector;
using std::set;
using std::array;

namespace jags {

static array<int, 2> mkDepth(vector<Node const *> const &parents)
{
    int sd = 0;
    for (unsigned long i = 0; i < parents.size(); ++i) {
	int sdi = parents[i]->depth()[0];
	if (sdi > sd) sd = sdi;
    }
    int dd = 1;
    for (unsigned long i = 0; i < parents.size(); ++i) {
	array<int, 2> const &pdepth = parents[i]->depth();
	if (pdepth[0] == sd && pdepth[1] >= dd) dd = pdepth[1] + 1;
    }
    
    array<int, 2> depth = {sd, dd};
    return depth;
}


DeterministicNode::DeterministicNode(vector<unsigned long> const &dim,
				     unsigned int nchain,
                                     vector<Node const *> const &parents)
    : Node(dim, nchain, parents), _fixed(true), _depth(mkDepth(parents))
{
    //Add this node as a deterministic child of its parents
    //taking care to avoid repeats
    set<Node const*> pset;
    for (unsigned long i = 0; i < parents.size(); ++i) {
	Node const *p = parents[i];
	if (pset.insert(p).second) {
	    p->addChild(this);
	}
    }

    //Deterministic nodes are not fixed if any parents are not fixed
    vector<Node const*>::const_iterator p;
    for (p = parents.begin(); p != parents.end(); ++p)
    {
	if (!(*p)->isFixed()) {
	    _fixed =  false;
	    break;
	}
    }

    /* 
       Fixed deterministic nodes should be immediately initialized by
       calling deterministicSample. We can't do it here because that
       is a virtual function. So we have to do it in whatever
       sub-class defines the deterministicSample member function.
    */
}

DeterministicNode::~DeterministicNode()
{
}

array<int, 2> const &DeterministicNode::depth() const {
   return _depth;
}
    
void DeterministicNode::randomSample(RNG*, unsigned int chain) {
    deterministicSample(chain);
}

RVStatus DeterministicNode::randomVariableStatus() const 
{ 
    return RV_FALSE; 
}

bool DeterministicNode::isFixed() const
{
    return _fixed;
}

    void DeterministicNode::unlinkParents()
    {
	for (unsigned long i = 0; i < parents().size(); ++i) {
	    parents()[i]->removeChild(this);
	}
    }
    
} //namespace jags
