#ifndef CONSTANT_FACTORY_H_
#define CONSTANT_FACTORY_H_

#include <graph/Graph.h>
#include <compiler/NodeFactory.h>

#include <map>
#include <cfloat>

class ConstantNode;

/**
 * Fuzzy "less than" operator for doubles, allowing some numerical
 * imprecision.
 */
bool lt(double arg1, double arg2);
/**
 * Function object for the STL map class using double as a key
 */
struct ltdouble
{
  bool operator()(double arg1, double arg2)
  {
    return lt(arg1, arg2);
  }
};

/**
 * @short Factory object for constant nodes 
 *
 * The purpose of a ConstantFactory is to avoid unnecessary
 * duplication of constant nodes by having a container class and
 * factory object for them that will create and/or lookup constant
 * nodes based on their value.
 */
class ConstantFactory : public NodeFactory
{ 
  unsigned int _nchain;
  std::map<double, ConstantNode*, ltdouble> _constmap;
 public:
  ConstantFactory(unsigned int nchain);
  /**
   * Get a constant node with a given value.  The results are cached,
   * so if a request is repeated, the same node will be returned.
   */
  ConstantNode *getConstantNode(double value);
};

#endif /* CONSTANT_FACTORY_H_ */
