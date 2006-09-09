#ifndef NODE_NAME_TAB_H_
#define NODE_NAME_TAB_H_

#include <string>

class Node;

/**
 * @short Lookup table for node names
 *
 * An abstract base class that provides a name for a given node.
 * @see Node##name
 */
class NodeNameTab {
 public:
  /**
   * Virtual destructor
   */
  virtual ~NodeNameTab();
  /**
   * Returns the name of the node, if it is included in the NodeNameTab,
   * or an empty string if not.
   */
  virtual std::string getName(Node const *node) const = 0;
};

#endif
