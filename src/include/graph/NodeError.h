#ifndef NODE_ERROR_H_
#define NODE_ERROR_H_

#include <stdexcept>

class Node;

/**
 * @short Exception class for Nodes
 */
class NodeError : public std::runtime_error {
 public:
  Node const * node;
  NodeError(Node const *enode, std::string const &msg);
};

#endif /* NODE_ERROR_H_ */
