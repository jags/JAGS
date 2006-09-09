#include <config.h>
#include <graph/MixtureNode.h>
#include <sarray/Index.h>

#include <utility>
#include <vector>
#include <cfloat>
#include <stdexcept>

#include <graph/NodeError.h>

using std::pair;
using std::vector;
using std::map;
using std::invalid_argument;
using std::logic_error;
using std::set;

static vector<Node*> mkParents(vector<Node *> const &index,
			       vector<pair<Index, Node*> > const &param)
{
  vector<Node*> parents;
  parents.reserve(index.size() + param.size());
  for (unsigned int i = 0; i < index.size(); ++i) {
    parents.push_back(index[i]);
  }
  for (unsigned int i = 0; i < param.size(); ++i) {
    parents.push_back(param[i].second);
  }
  return parents;
}

#include <iostream>
MixtureNode::MixtureNode(vector<Node *> const &index,
			 vector<pair<Index, Node *> > const &parameters)
  : DeterministicNode(parameters[0].second->dim(true),
		      mkParents(index, parameters)),
		      _Nindex(index.size())
{
  for (vector<Node *>::const_iterator i = index.begin(); i != index.end(); ++i)
    {
      Node *node = *i;
      if (node->length() != 1 || !node->isDiscreteValued()) {
         std::cout << "length = " << node->length() << "\n";
         throw NodeError(node, "Invalid index parameter for mixture node");
         //throw invalid_argument("Invalid index parameter for mixture node");
      }
    }

  unsigned int ndim = parameters.size();
  Index const &default_dim = dim(false);
  for (unsigned int i = 0; i < ndim; ++i) {
    Node *node = parameters[i].second;
    if (!node) {
      throw invalid_argument("Null parameter in MixtureNode");
    }
    _map[parameters[i].first] = node;
  }
  
  bool isdiscrete = true;
  for (unsigned int i = 0; i < ndim; ++i) {
    Node *node = parameters[i].second;
    if (node->dim(true) != default_dim) {
      throw invalid_argument("Dimension mismatch for MixtureNode parameters");
    }
    if (!node->isDiscreteValued()) {
      isdiscrete = false;
    }
  }
  if (isdiscrete) {
    setDiscreteValued();
  }
}

MixtureNode::~MixtureNode()
{
}

//debuggin
#include <iostream>
#include <sarray/Range.h>
#include <sarray/SArray.h>
void MixtureNode::deterministicSample(unsigned int chain)
{
  Index i(_Nindex);
  vector <Node*> const &parents = this->parents();
  for (unsigned int j = 0; j < _Nindex; ++j) {
    i[j] = static_cast<long>(parents[j]->value(chain)[0]);
  }
  map<Index,Node*>::iterator p = _map.find(i);
  if (p != _map.end()) {
    setValue(p->second->value(chain), p->second->length(), chain);
  }
  else {
    std::cout << "Got " << print(i) << "\nOriginally\n";
    for (unsigned int j = 0; j < _Nindex; ++j) {
        std:: cout << parents[j]->value(chain)[0] << "\n";
        if (parents[j]->value(chain)[0] == JAGS_NA)
           std::cout << "(which is  missing)\n";
    }
    std::cout << "Expected one of \n";
    for (p = _map.begin(); p != _map.end(); ++p) {
       std::cout << print(p->first) << "\n";
    }
    throw logic_error("Invalid index in MixtureNode");
  }
}

unsigned int MixtureNode::index_size() const
{
  return _Nindex;
}

MixtureNode const *asMixture(Node const *node)
{
  return dynamic_cast<MixtureNode const*>(node);
}

bool isMixture(Node const *node)
{
  return dynamic_cast<MixtureNode const*>(node);
}

bool MixtureNode::isLinear(set<Node*> const &parameters, bool fixed) const
{
  if (fixed)
      return false;

  //Check that none of the indices are in the parameter set
  vector<Node*> const &par = parents();
  for (unsigned int i = 0; i < _Nindex; ++i) {
    if (parameters.count(par[i]))
      return false;
  }
  return true;
}

bool MixtureNode::isScale(set<Node*> const &parameters, bool fixed) const
{
    return isLinear(parameters, fixed);
}
