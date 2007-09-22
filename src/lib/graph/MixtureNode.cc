#include <config.h>
#include <graph/MixtureNode.h>
#include <graph/NodeNameTab.h>

#include <utility>
#include <vector>
#include <stdexcept>

#include <graph/NodeError.h>

using std::vector;
using std::map;
using std::invalid_argument;
using std::logic_error;
using std::set;
using std::string;

/*
  Calculates the dimensions of a mixture node given its possible parent
  values. If the parents have inconsistent dimensions, then a logic
  error is thrown.
*/
static vector<unsigned int> const &
mkDim(map<vector<int>, Node const *> const &mixmap)
{
    map<vector<int>, Node const *>::const_iterator p = mixmap.begin();
    vector<unsigned int> const &dim = p->second->dim();
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
mkParents(vector<Node const *> const &index,
	  map<vector<int>, Node const *> const &mixmap)
{
    vector<Node const *> parents;
    parents.reserve(index.size() + mixmap.size());
    for (unsigned int i = 0; i < index.size(); ++i) {
	parents.push_back(index[i]);
    }
    for (map<vector<int>, Node const *>::const_iterator p = mixmap.begin();
	 p != mixmap.end(); ++p) 
    {
	parents.push_back(p->second);
    }
    return parents;
}

MixtureNode::MixtureNode (vector<Node const *> const &index,
			  map<vector<int>, Node const *> const &mixmap)
    : DeterministicNode(mkDim(mixmap), mkParents(index, mixmap)),
      _map(mixmap), _Nindex(index.size())
{
    // Check validity of index argument

    if (index.empty())
	throw invalid_argument("NULL index in MixtureNode constructor");

    for (vector<Node const *>::const_iterator i = index.begin(); 
	 i != index.end(); ++i)
    {
	Node const *node = *i;
	if (node->length() != 1 || !node->isDiscreteValued() ||
	    node->isObserved()) 
	{
	    throw invalid_argument("Invalid index in MixtureNode constructor");
	}
    }

    // Check validity of mixmap argument

    if (mixmap.size() < 2)
	throw invalid_argument("Trivial mixmap in MixtureNode constructor");

    for (map<vector<int>, Node const *>::const_iterator p = mixmap.begin();
	 p != mixmap.end(); ++p)
    {
	if (p->first.size() != _Nindex) {
	    throw invalid_argument("Invalid index in MixtureNode");
	}
    }

    // A mixture node is discrete if all its parents are discrete

    bool isdiscrete = true;
    vector<Node const *> const &par = parents();
    for(unsigned int i = _Nindex; i < par.size(); ++i) {
	if (!par[i]->isDiscreteValued()) {
	    isdiscrete = false;
	    break;
	}
    }

    if (isdiscrete) {
	setDiscreteValued();
    }
}

MixtureNode::~MixtureNode()
{
}

/* Do not delete commented sections: they are useful for debugging
#include <iostream>
#include <sarray/Range.h>
#include <sarray/nainf.h>
#include <graph/NodeError.h>
*/
void MixtureNode::deterministicSample(unsigned int chain)
{
    vector<int> i(_Nindex);
    vector <Node const*> const &par = parents();
    for (unsigned int j = 0; j < _Nindex; ++j) {
	i[j] = static_cast<int>(*par[j]->value(chain));
    }
    map<vector<int>, Node const *>::const_iterator p = _map.find(i);
    if (p == _map.end()) {
	/*
	std::cout << "Got " << print(Range(i)) << "\nOriginally\n";
	for (unsigned int j = 0; j < _Nindex; ++j) {
	    std::cout << par[j]->value(chain)[0] << "\n";
	    if (par[j]->value(chain)[0] == JAGS_NA)
		std::cout << "(which is  missing)\n";
	}
	std::cout << "Expected one of \n";
	for (p = _map.begin(); p != _map.end(); ++p) {
	    std::cout << print(Range(p->first)) << "\n";
	}
	*/
	throw NodeError(this, "Invalid index in MixtureNode");
    }
    else {
	setValue(p->second->value(chain), length(), chain);	
    }
    }

unsigned int MixtureNode::index_size() const
{
  return _Nindex;
}

string MixtureNode::name(NodeNameTab const &name_table) const
{
    string name = name_table.getName(this);
    if (!name.empty())
	return name;

    name = "mixture(index=[";

    vector<Node const *> const &par = parents();
    vector<int> i(_Nindex);
    for (unsigned int j = 0; j < _Nindex; ++j) {
	if (j > 0) {
	    name.append(",");
	}
	name.append(par[j]->name(name_table));
    }
    name.append("], parents= ");

    /* We can't list all possible parents in a name, since there is
       no limit on the number. So we take the first and last */
    name.append(par[_Nindex]->name(name_table)); //first parent
    if (par.size() > _Nindex + 2) {
        name.append("...");
    }
    else {
        name.append(",");
    }
    name.append(par[par.size() - 1]->name(name_table)); //last parent
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

bool 
MixtureNode::isLinear(set<Node const*> const &parameters, bool fixed) const
{
  if (fixed)
      return false;

  //Check that none of the indices are in the parameter set
  vector<Node const*> const &par = parents();
  for (unsigned int i = 0; i < _Nindex; ++i) {
    if (parameters.count(par[i]))
      return false;
  }
  return true;
}

bool MixtureNode::isScale(set<Node const*> const &parameters, bool fixed) const
{
    return isLinear(parameters, fixed);
}

bool MixtureNode::checkParentValues(unsigned int chain) const
{
    return true;
}
