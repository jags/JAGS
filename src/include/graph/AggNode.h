#ifndef AGGREGATE_NODE_H_
#define AGGREGATE_NODE_H_

#include <graph/DeterministicNode.h>
#include <vector>

/**
 * An aggregate Node copies its values directly from its parents.  As
 * the name suggests, an aggregate node can combine several smaller
 * nodes, but in fact the Aggregate node is very general and can take
 * subsets of the values of its parents.
 *
 * @short Aggregate Node
 */
class AggNode : public DeterministicNode {
  std::vector<unsigned long> _offsets;

  /* Forbid copying */
  AggNode(AggNode const &orig);
  AggNode &operator=(AggNode const &rhs);
 public:
  /**
   * Constructor.  
   *
   * @param dim Dimension of the Node.
   *
   * @param nodes Vector of parent Nodes of size dim.size(). A node
   * may appear several times in the vector.
   *
   * @param offsets Vector of offsets of size dim.length().  The
   * offset denotes the element of the value vector to be copied
   * from the corresponding parent.
   */
  AggNode(Index const &dim,
	  std::vector<Node*> const &nodes, 
	  std::vector<unsigned long> const &offsets);
  ~AggNode();

  /**
   * Copies values from parents.
   */
  void deterministicSample(unsigned int chain);
  /**
   * An AggNode always preserves linearity
   */
  bool isLinear(std::set<Node*> const &parameters, bool fixed) const;
  /**
   * An AggNode is a scale transformation only if it has a single
   * parent.
   */
  bool isScale(std::set<Node*> const &parameters, bool fixed) const;
};

AggNode const *asAggregate(Node *node);

#endif /* AGGREGATE_NODE_H */
