#ifndef LOGICAL_NODE_H_
#define LOGICAL_NODE_H_

#include <graph/DeterministicNode.h>
#include <vector>

namespace jags {

class Function;

/**
 * @short Node defined by the BUGS-language operator <-
 *
 * A LogicalNode is a type of DeterministicNode that is defined by a
 * function, and a vector of parents as parameters. In the BUGS
 * language the function may be represented as an operator.
 */
class LogicalNode : public DeterministicNode {
    Function const * const _func;
    bool _discrete;
protected:
    const std::vector<std::vector<double const*> > _parameters;
    void initializeFixed();
public:
    /**
     * A logical node is defined by a function (which may be an inline
     * operator in the BUGS language) and its parameters.
     */
    LogicalNode(std::vector<unsigned int> const &dim,
		unsigned int nchain,
                std::vector<Node const*> const &parameters,
		Function const *func);
    /**
     * A LogicalNode may be discrete valued if its parents are. 
     * @see Function##isDiscreteValued
     */
    bool isDiscreteValued() const;
    /**
     * A LogicalNode belongs to a closed class if its function
     * preserves that class.
     *
     * @see Function#isLinear, Function#isScale, Function#isPower
     */
    bool isClosed(std::set<Node const *> const &ancestors, 
		  ClosedFuncClass fc, bool fixed) const;
    std::string deparse(std::vector<std::string> const &) const;
};

} /* namespace jags */

#endif /* LOGICAL_NODE_H_ */
