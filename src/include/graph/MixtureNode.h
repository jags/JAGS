#ifndef MIXTURE_NODE_H_
#define MIXTURE_NODE_H_

#include <graph/DeterministicNode.h>
#include <map>
#include <vector>

/**
 * A mixture node takes its value from one of several parents, based
 * on the current value of an index node, which is also a parent
 * @short Node for mixture models
 */
class MixtureNode : public DeterministicNode {
  std::map<Index, Node*> _map;
  unsigned int _Nindex;
public:
  /**
   * Constructor. 
   * @param index Vector of discrete-valued scalar nodes
   * @param parameters Vector of pairs. Each pair associates a
   * possible value of the index nodes with a parent. The mixture node
   * copies its value from the parent that matches the current value
   * of the index.
   */
  MixtureNode(std::vector<Node *> const &index,
	      std::vector<std::pair<Index, Node *> > const &parameters);
  ~MixtureNode();
	
  void deterministicSample(unsigned int chain);
  /** 
   * Returns the vector of index nodes.
   */
  //std::vector<Node *> const &index() const;
  /**
   * Returns true if all of the parameters of the node are discrete-valued
   */
  //FIXME: This isn't a virtual function
  //bool isDiscreteValued() const;
  /**
   * Returns the number of index parameters
   */
  unsigned int index_size() const;
  /**
   * A MixtureNode preserves linearity if none of the indices are
   * parameters. It is never a fixed linear function.
   */
  bool isLinear(std::set<Node*> const &parameters, bool fixed) const;
  /**
   * A MixtureNode is a scale transformation if none of the indices
   * are parameters. It is never a fixed linear function.
   */
  bool isScale(std::set<Node*> const &parameters, bool fixed) const;
};

MixtureNode const *asMixture(Node const *node);
bool isMixture(Node const *node);

#endif /* MIXTURE_NODE_H_ */
