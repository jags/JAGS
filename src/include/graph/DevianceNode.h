#ifndef DEVIANCE_NODE_H_
#define DEVIANCE_NODE_H_

#include <graph/DeterministicNode.h>
#include <graph/StochasticNode.h>
#include <vector>
#include <set>
#include <stdexcept>

class Model;
class StochasticNode;

class DevianceNode : public DeterministicNode
{
  std::vector<StochasticNode const*> _parameters;
 public:
  DevianceNode(std::set<StochasticNode*> const &parameters);
  void deterministicSample(unsigned int chain);
  /**
   * A deviance node never preserves linearity
   */
  bool isLinear(std::set<Node*> const &parameters, bool fixed) const;
  /**
   * A deviance node is never a scale transformation of its parameters
   */
  bool isScale(std::set<Node*> const &parameters, bool fixed) const;
};


#endif
