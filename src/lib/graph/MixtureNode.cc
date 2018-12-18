#include <config.h>
#include <graph/MixtureNode.h>
#include <graph/GraphMarks.h>
#include <graph/Graph.h>
#include <graph/MixTab.h>
#include <graph/NodeError.h>

//Debugging only
//#include <util/nainf.h>

#include <utility>
#include <vector>
#include <stdexcept>

using std::vector;
using std::map;
using std::invalid_argument;
using std::logic_error;
using std::set;
using std::string;
using std::pair;


namespace jags {

    typedef map<MixMap, pair<MixTab*, int> > MixTabMap;

    static MixTabMap &mixTabMap()
    {
	// Repository of MixTab objects that are shared between
	// MixtureNodes.  We use reference counting to keep track of
	// whether a MixTab is in use.
	static MixTabMap _map;
	return _map;
    }

    static MixTab const *getTable(MixMap const &mixmap)
    {
	// Returns a MixTab object from the repository corresponding
	// to the given MixMap.

	// N.B. This must be called only by the MixtureNode
	// constructor!

	MixTabMap &tabmap = mixTabMap();
	MixTabMap::iterator p = tabmap.find(mixmap);
	if (p == tabmap.end()) {
	    //MixTab does not exist in repository; create and
	    //insert a new one.
	    MixTab *newtab = new MixTab(mixmap);
	    MixTabMap::mapped_type newentry(newtab, 1);
	    p = tabmap.insert(MixTabMap::value_type(mixmap, newentry)).first;
	}
	else {
	    //MixTab already exists in the repository; increment
	    //reference count
	    p->second.second++;
	}
	return p->second.first;
    }

    static MixTabMap::iterator findTable(MixTab const *table)
    {
	// Inverse lookup of the MixTab repository using the MixTab

	MixTabMap &tabmap = mixTabMap();

	MixTabMap::iterator p = tabmap.begin();
	
	for( ; p != tabmap.end(); ++p) {
	    if (p->second.first == table) {
		break;
	    }
	}
	if (p == tabmap.end()) {
	    throw logic_error("Failed to find MixTab in MixtureNode");
	}
	return p; //Wall

    }

    static void removeTable(MixTab const *table)
    {
	// Decrements the reference count of a MixTab in the shared
	// repository.  When the reference count reaches zero, the MixTab
	// is removed.
	
	// N.B. This must be called only by the MixtureNode destructor!
	
	MixTabMap::iterator p = findTable(table);
	p->second.second--;
	if (p->second.second == 0) {
	    mixTabMap().erase(p);
	}
    }
    
    /*
      Calculates the dimensions of a mixture node given its possible
      parent values. If the parents have inconsistent dimensions, then
      a logic error is thrown.
    */
    static vector<unsigned long> const &mkDim(MixMap const &mixmap)
    {
	MixMap::const_iterator p = mixmap.begin();
	vector<unsigned long> const &dim = p->second->dim();
	for (++p ; p != mixmap.end(); ++p) {
	    if (p->second->dim() != dim) {
		throw logic_error("Dimension mismatch in MixtureNode parents");
	    }
	}
	return dim;
    }
    
    /* 
       Creates a vector of parent nodes from the arguments passed to the
       constructor. 
       
       The index nodes come first, in the order supplied, then the parents
       supplied in the mixmap parameter, in the order determined by the
       corresponding indices.
    */
    static vector<Node const *> 
    mkParents(vector<Node const *> const &index, MixMap const &mixmap)
    {
	vector<Node const *> parents;
	parents.reserve(index.size() + mixmap.size());
	for (unsigned long i = 0; i < index.size(); ++i) {
	    parents.push_back(index[i]);
	}
	for (map<vector<unsigned long>, Node const *>::const_iterator p = mixmap.begin();
	     p != mixmap.end(); ++p) 
	{
	    parents.push_back(p->second);
	}
	return parents;
    }

MixtureNode::MixtureNode (vector<Node const *> const &index,
			  unsigned int nchain,
			  MixMap const &mixmap)
    : DeterministicNode(mkDim(mixmap), nchain, mkParents(index, mixmap)),
      _table(getTable(mixmap)), _Nindex(index.size()), _discrete(true),
      _active_parents(nchain)
{
    // Check validity of index argument

    if (index.empty())
	throw invalid_argument("NULL index in MixtureNode constructor");

    for (vector<Node const *>::const_iterator i = index.begin(); 
	 i != index.end(); ++i)
    {
	Node const *node = *i;
	if (node->length() != 1 || !node->isDiscreteValued() ||
	    node->isFixed()) 
	{
	    throw invalid_argument("Invalid index in MixtureNode constructor");
	}
    }

    // Check validity of MixMap argument
    if (mixmap.size() < 2)
	throw invalid_argument("Trivial MixMap in MixtureNode constructor");

    // Check consistency of arguments
    if (index.size() != mixmap.begin()->first.size()) {
	throw invalid_argument("Dimension mismatch in MixtureNode constructor");
    }

    //Check discreteness 
    vector<Node const *> const &par = parents();
    for (unsigned long i = _Nindex; i < par.size(); ++i)
    {
	if (!par[i]->isDiscreteValued()) {
	    _discrete = false;
	    break;
	}
    }
}

MixtureNode::~MixtureNode()
{
    removeTable(_table);
}

/* Do not delete commented sections: they are useful for debugging
#include <iostream>
#include <sarray/Range.h>
#include <util/nainf.h>
#include <graph/NodeError.h>
*/

void MixtureNode::updateActive(unsigned int chain)
{
    vector<unsigned long> i(_Nindex);
    vector <Node const*> const &par = parents();
    for (unsigned long j = 0; j < _Nindex; ++j) {
	i[j] = static_cast<unsigned long>(*par[j]->value(chain));
    }

    _active_parents[chain] = _table->getNode(i);
    if (_active_parents[chain] == nullptr) {
	/*
	std::cout << "Got " << printIndex(i) << "\nOriginally\n";
	for (unsigned int j = 0; j < _Nindex; ++j) {
	    std::cout << par[j]->value(chain)[0] << "\n";
	    if (par[j]->value(chain)[0] == JAGS_NA)
		std::cout << "(which is  missing)\n";
	}
	*/
	string msg = string("Invalid index") + printIndex(i) +
	    " in mixture node";
	throw NodeError(this, msg);
    }
}

void MixtureNode::deterministicSample(unsigned int chain)
{
    updateActive(chain);
    setValue(_active_parents[chain]->value(chain), length(), chain);	
}

Node const *MixtureNode::activeParent(unsigned int chain) const
{
    return _active_parents[chain];
}

unsigned long MixtureNode::index_size() const
{
  return _Nindex;
}

string MixtureNode::deparse(vector<string> const &parents) const
{
    string name = "mixture(index=[";

    vector<unsigned long> i(_Nindex);
    for (unsigned long j = 0; j < _Nindex; ++j) {
	if (j > 0) {
	    name.append(",");
	}
	name.append(parents[j]);
    }
    name.append("], parents= ");

    /* We can't list all possible parents in a name, since there is
       no limit on the number. So we take the first and last */
    name.append(parents[_Nindex]); //first parent
    if (parents.size() > _Nindex + 2) {
        name.append("...");
    }
    else {
        name.append(",");
    }
    name.append(parents.back()); //last parent
    name.append(")");
	      
    return name;
}

MixtureNode const *asMixture(Node const *node)
{
  return dynamic_cast<MixtureNode const*>(node);
}

bool isMixture(Node const *node)
{
  return dynamic_cast<MixtureNode const*>(node);
}

bool MixtureNode::isDifferentiable(Node const *arg) const
{
    auto par = parents();
    for (unsigned long i = _Nindex; i < par.size(); ++i) {
	if (par[i]==arg) return false;
    }
    return true;
}

void MixtureNode::gradient(double *grad, Node const *arg,
			   unsigned int chain) const
{
    if (arg == _active_parents[chain]) {
	for (unsigned int i = 0; i < _length; ++i) {
	    grad += 1;
	}
    }
}

bool MixtureNode::isClosed(set<Node const *> const &ancestors, 
			   ClosedFuncClass fc, bool fixed) const
{
    if (fixed)
	return false;

    //Check that none of the indices are in the ancestor set
    vector<Node const*> const &par = parents();
    for (unsigned long i = 0; i < _Nindex; ++i) {
	if (ancestors.count(par[i])) {
	    return false;
	}
    }

    switch(fc) {
    case DNODE_LINEAR: case DNODE_SCALE_MIX: case DNODE_POWER:
	break;
    case DNODE_SCALE: case DNODE_ADDITIVE:
	//Only a scale or additive function if all possible parents are scale
	//or additive functions, respectively.
	for (unsigned long i = _Nindex; i < par.size(); ++i) {
	    if (ancestors.count(par[i])==0)
		return false;
	}
	break;
    }
    
    return true;
}

bool MixtureNode::checkParentValues(unsigned int) const
{
    return true;
}

    /*
DeterministicNode *MixtureNode::clone(vector<Node const *> const &parents) const
{
    vector<Node const *> index(_Nindex);
    vector<Node const *>::const_iterator p = parents.begin();
    for (unsigned int i = 0; i < _Nindex; ++i) {
	index[i] = *p;
	++p;
    }
    
    // Find the MixMap corresponding to this node, and copy it
    MixMap mixmap = findTable(_table)->first;
    // Replace entries in the copy
    MixMap::const_iterator q = mixmap.begin();
    while (p != parents.end() && q != mixmap.end()) {
	mixmap[q->first] = *p;
	++q;
	++p;
    }
    
    return new MixtureNode(index, mixmap);
}
    */
    
bool MixtureNode::isDiscreteValued() const
{
    return _discrete;
}

    MixTab const *MixtureNode::mixTab() const
    {
	return _table;
    }

} //namespace jags
