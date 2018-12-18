#ifndef ARRAY_LOGICAL_NODE_H_
#define ARRAY_LOGICAL_NODE_H_

#include <graph/LogicalNode.h>
#include <vector>

namespace jags {

class ArrayFunction;

/**
 * @short Array-valued Node defined by the BUGS-language operator <-
 */
class ArrayLogicalNode : public LogicalNode {
    ArrayFunction const * const _func;
    std::vector<std::vector<unsigned long> >  _dims;
public:
    /**
     * A logical node is defined by a function (which may be an inline
     * operator in the BUGS language) and its parameters.
     */
    ArrayLogicalNode(ArrayFunction const *func, unsigned int nchain,
		     std::vector<Node const*> const &parameters);
    /**
     * Calculates the value of the node based on the parameters. 
     */
    void deterministicSample(unsigned int chain);
    /**
     * @see ArrayFunction#checkParameterValue.
     */
    bool checkParentValues(unsigned int chain) const;
    //DeterministicNode *clone(std::vector<Node const *> const &parents) const;
    void gradient(double *grad, Node const *arg, unsigned int chain) const;
};

} /* namespace jags */

#endif /* ARRAY_LOGICAL_NODE_H_ */
