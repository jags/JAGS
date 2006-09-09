#include <config.h>
#include <compiler/MixtureFactory.h>
#include <graph/MixtureNode.h>
#include <model/NodeArray.h>

#include <stdexcept>
#include <utility>

using std::pair;
using std::map;
using std::vector;
using std::runtime_error;
using std::logic_error;

bool compMixPair(MixPair const &arg1, MixPair const &arg2)
{
  //First compare indices
  vector<Node *> const &index1 = arg1.first;
  vector<Node *> const &index2 = arg2.first;
  if (index1.size() < index2.size()) {
    return true;
  }
  else if (index1.size() > index2.size()) {
    return false;
  }
  else {
    for (unsigned int i = 0; i < index1.size(); ++i) {
      if (index1[i] < index2[i]) {
	return true;
      }
      else if (index1[i] > index2[i]) {
	return false;
      }
    }
  }

  //Same indices. Now compare parameters 
  vector<pair<Index, Node*> > const &parameters1 = arg1.second;
  vector<pair<Index, Node*> > const &parameters2 = arg2.second;
  if (parameters1.size() < parameters2.size()) {
    return true;
  }
  else if (parameters1.size() > parameters2.size()) {
    return false;
  }
  else {
    for (unsigned int i = 0; i < parameters1.size(); ++i) {
      if (parameters1[i].first < parameters2[i].first) {
	return true;
      }
      else if (parameters2[i].first < parameters1[i].first) {
	return false;
      }
      else if (parameters1[i].second < parameters2[i].second) {
	return true;
      }
      else if (parameters1[i].second > parameters2[i].second) {
	return false;
      }
    }
    return false; //equal
  }
}

MixtureNode *
MixtureFactory::getMixtureNode(vector<Node*> const &index,
			       vector<pair<Index, Node*> > const &parameters)
{
  MixPair mpair(index, parameters);
  map<MixPair, MixtureNode*, ltmixpair>::const_iterator p = 
    _mixmap.find(mpair);

  if (p != _mixmap.end()) {
    return p->second;
  }
  else {
    MixtureNode *mix = new MixtureNode(index, parameters);
    _mixmap[mpair] = mix;
    _graph.add(mix);
    return mix;
  }
}
