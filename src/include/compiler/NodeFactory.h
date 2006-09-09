#ifndef NODE_FACTORY_H_
#define NODE_FACTORY_H_

#include <graph/Graph.h>

/**
 * @short Abstract factory class for deterministic nodes The value of
 * any deterministic node is determined uniquely by its parents.  The
 * purpose of a NodeFactory is to avoid unnecessary duplication of
 * deterministic nodes by providing a container class and factory
 * object for them that will create and/or lookup deterministic nodes.
 */
class NodeFactory { 
 protected:
  Graph _graph;
 public:
  /**
   * Returns the graph of nodes created by the NodeFactory
   */
  virtual ~NodeFactory();
  Graph const &graph() const;
};

/**
 * Two doubles are considered equal if their values lie within
 * a certain numerical tolerance.
 */
bool equal(double arg1, double arg2);
/**
 * Two nodes are considered equal if they have the same address, or if
 * they are both constant and have the same value (within numerical
 * tolerance), or if they are deterministic, fixed nodes with the same
 * value.
*/
bool equal(Node const* node1, Node const* node2);

#endif /* NODE_FACTORY_H_ */
