#ifndef DETERMINISTIC_NODE_H_
#define DETERMINISTIC_NODE_H_

#include <graph/Node.h>

class DeterministicNode : public Node {
 public:
  DeterministicNode(Index const &dim, std::vector<Node*> const &parents);
  void randomSample(RNG*, unsigned int nchain) {deterministicSample(nchain);};
  /**
   * Deterministic nodes are not random variables
   */
  bool isVariable() const {return false;}
  /**
   * Tests whether the value of the node is a linear function of the
   * ancestor nodes X1, X2, ... Xn of form a + b %*% X1 + b2 %*% X2 +
   * ...  + bn %*% Xn.
   *
   * False negative responses are permitted: i.e. the value false may
   * be returned when the node is, in fact, a linear function, but
   * false positives are not allowed.
   * 
   * @param parameters Set of all nodes that are themselves linear
   * functions functions of the ancestor nodes X1, X2, ... Xn.  The
   * set must include the ancestor nodes (as they are trivial linear
   * functions of themselves).  The current node may be a non-linear
   * function of nodes not in this set.
   *
   * @param fixed When true, the test is more stringent and returns
   * the value true is returned only if the function is linear and the
   * coefficients b1, ... bn are fixed (but not necessarily a).
   */
  virtual bool 
    isLinear(std::set<Node*> const &parameters, bool fixed) const = 0;
  /**
   * Tests whether the value of the node is a function of ancestor
   * node X of the form b %*% X. False negative responses are
   * permitted.
   * 
   * @param parameters Set of all nodes that are themselves scale
   * transformations of the ancestor node X. The set must include X
   * (as this is a trivial scale transformation of itself).
   */
  virtual bool 
      isScale(std::set<Node*> const &parameters, bool fixed) const = 0;
};

#endif /* DETERMINISTIC_NODE_H_ */
