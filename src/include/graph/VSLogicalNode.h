#ifndef V_SCALAR_LOGICAL_NODE_H_
#define V_SCALAR_LOGICAL_NODE_H_

#include <graph/LogicalNode.h>
#include <vector>

namespace jags {

class ScalarFunction;

/**
 * @short Vectorized node using a scalar function
 *
 * A VSLogicalNode is a logical node that uses a ScalarFunction as
 * the defining function. The node is vectorized.
 */
class VSLogicalNode : public LogicalNode {
    ScalarFunction const * const _func;
    std::vector<bool> _isvector;
public:
    /**
     * A scalar logical node is defined by a scalar function (which
     * may be an inline operator in the BUGS language) and its
     * parameters.
     */
    VSLogicalNode(ScalarFunction const *func, unsigned int nchain,
		  std::vector<Node const*> const &parameters);
    /**
     * Calculates the value of the node based on the parameters. 
     */
    void deterministicSample(unsigned int chain);
    /**
     * @see ScalarFunction#checkParameterValue.
     */
    bool checkParentValues(unsigned int chain) const;
    //DeterministicNode *clone(std::vector<Node const *> const &parents) const;
    void gradient(double *grad, Node const *arg, unsigned int chain) const;
};

} /* namespace jags */

#endif /* V_SCALAR_LOGICAL_NODE_H_ */
