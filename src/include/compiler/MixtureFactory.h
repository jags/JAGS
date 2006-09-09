#ifndef MIXTURE_FACTORY_H_
#define MIXTURE_FACTORY_H_

#include <graph/Graph.h>
#include <graph/Node.h>
#include <sarray/SArray.h>
#include <compiler/NodeFactory.h>

#include <vector>
#include <map>
#include <cfloat>

class NodeArray;
class MixtureNode;

typedef std::pair<std::vector <Node*>, std::vector<std::pair<Index, Node*> > > 
MixPair;

bool compMixPair(MixPair const &, MixPair const &);

struct ltmixpair
{
  bool operator()(MixPair const &arg1, MixPair const &arg2)
  {
    return compMixPair(arg1, arg2);
  }
};

/**
 * @short Factory object for mixture nodes 
 * 
 * The purpose of a MixtureFactory is to avoid unnecessary duplication
 * of mixture nodes by having a container class and factory object
 * that will create and/or lookup mixture nodes.
 */
class MixtureFactory : public NodeFactory
{ 
  std::map<MixPair, MixtureNode*, ltmixpair> _mixmap;
public:
  /**
   * Get a mixture node.  The results are cached, so if a request is
   * repeated, the same node will be returned.
   *
   * @param index Vector of discrete-valued Nodes that index the
   * possible parameters
   * @param parameters Vector of pairs. The Index shows a possible value
   * of the index vector and the Node pointer shows the corresponding
   * parent from which the mixture node copies its value.
   *
   */
  MixtureNode *getMixtureNode(std::vector<Node *> const &index,
			      std::vector<std::pair<Index, Node*> > 
			      const &parameters);
};

#endif /* MIXTURE_FACTORY_H_ */
