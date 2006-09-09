#ifndef CONSTANT_NODE_H_
#define CONSTANT_NODE_H_

#include <graph/Node.h>

struct RNG;

/**
 * Constant nodes represent the top level parameters in any directed
 * graphical model.  They are always scalar, and should not have
 * parents.
 */
class ConstantNode : public Node {
 public:
  /**
   * Construct node and set its value. The value is fixed and is
   * shared between all chains.
   */
  ConstantNode(double value, unsigned int nchain);
  /**
   * This function does nothing. It exists only so that objects of
   * class ConstantNode can be instantiated.
   */
  void deterministicSample(unsigned int) {};
  void randomSample(RNG*, unsigned int) {};
  std::string name(NodeNameTab const &name_table) const;
  /**
   * Constant nodes are observed random variables whose values we are
   * conditioning on.
   */
  bool isVariable() const { return true;}
};

#endif /* CONSTANT_NODE_H_ */




