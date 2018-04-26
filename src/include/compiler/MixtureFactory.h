#ifndef MIXTURE_FACTORY_H_
#define MIXTURE_FACTORY_H_

#include <vector>
#include <map>
#include <cfloat>

#include <graph/MixtureNode.h>

namespace jags {

class NodeArray;
class Model;

/**
 * A "mixture pair" uniquely indexes a mixture node. The first element
 * is a vector of (discrete-valued) index nodes. The second element is
 * a vector of parent nodes from which the mixture node copies its
 * value.
 *
 * It is implicit in the mixture pair definition that the index nodes
 * take a range of possible values. If these values are calculated and
 * sorted, then ther there is a one-to-one correspondence between the
 * sorted index values and the vector of parent nodes.
 *
 * The reason this vector of possible index values is not used is
 * because it is not necessary to uniquely define a mixture node.
 */
typedef std::pair<std::vector<Node const*>, std::vector<Node const*> > MixPair;

/**
 * @short Factory for MixtureNode objects
 * 
 * The purpose of a MixtureFactory is to avoid unnecessary duplication
 * of mixture nodes by having a container class and factory object
 * that will create and/or lookup mixture nodes.
 */
class MixtureFactory  
{ 
  std::map<MixPair, MixtureNode*> _mix_node_map;
public:
  /**
   * Get a mixture node.  The results are cached, so if a request is
   * repeated, the same node will be returned.
   *
   * @param index Vector of discrete-valued Nodes that index the
   * possible parameters 
   *
   * @param parameters Vector of pairs. The first element of the pair
   * shows a possible value of the index vector and the second element
   * shows the corresponding parent from which the mixture node copies
   * its value.
   *
   * @param model Model to which newly allocated mixturenodes are added.
   *
   */
  MixtureNode *getMixtureNode(std::vector<Node const *> const &index,
			      MixMap const &parameters, Model &model);
};

} /* namespace jags */

#endif /* MIXTURE_FACTORY_H_ */
