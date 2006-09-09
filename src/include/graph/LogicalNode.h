#ifndef LOGICAL_NODE_H_
#define LOGICAL_NODE_H_

#include <graph/DeterministicNode.h>
#include <vector>

class Function;
class SArray;

/**
 * The value of a logical node is a function of its parents.
 *
 * @short Logical Node
 */
class LogicalNode : public DeterministicNode {
  Function const * const _func;
  std::vector<std::vector<SArray const*> > _parameters;

  /* Forbid copying of logical nodes */
  LogicalNode(LogicalNode const &orig);
  LogicalNode &operator=(LogicalNode const &rhs);

 public:
  LogicalNode(Function const *func, std::vector<Node*> const &parameters);
  ~LogicalNode();
  /** Returns the parameters of the node as a vector of SArray pointers  */
  std::vector<SArray const *> const &parameters(unsigned int chain) const;
  /** Returns the Function of the logical node */
  Function const* function() const;
  void deterministicSample(unsigned int chain);
  /** 
   * Calculates the value of the node based on the parameters. This
   * function calls Function::checkParameterValue and will throw an
   * exception if the result is false.
   */
  std::string name(NodeNameTab const &name_table) const;
  bool isLinear(std::set<Node*> const &parameters, bool fixed) const;
  bool isScale(std::set<Node*> const &parameters, bool fixed) const;
};

bool isLogical(Node const *node);
LogicalNode const *asLogical(Node const *node);

#endif /* LOGICAL_NODE_H_ */
