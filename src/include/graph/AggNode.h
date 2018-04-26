#ifndef AGGREGATE_NODE_H_
#define AGGREGATE_NODE_H_

#include <graph/DeterministicNode.h>
#include <util/logical.h>

namespace jags {

/**
 * @short Aggregate Node combining values from other nodes
 *
 * An aggregate Node copies its values directly from its parents: it can
 * be used to aggregate several small nodes into a larger one, or to take
 * a subset of a larger node, or some combination of the two.
 *
 */
class AggNode : public DeterministicNode {
    std::vector<unsigned long> _offsets;
    std::vector<double const *>  _parent_values;
    bool _discrete;
    /* Forbid copying */
    AggNode(AggNode const &orig);
    AggNode &operator=(AggNode const &rhs);
public:
    /**
     * The value vector of an AggNode satisfies the equality.
     * <pre>
     * value(chain)[i] == parents[i]->value(chain)[offsets[i]]
     * </pre>
     *
     * @param dim Dimension of the Node.
     *
     * @param nchain Number of chains
     * 
     * @param parents Vector of parent Nodes. This vector should have
     * the same size as the value array of the AggNode to be
     * constructed (or a length_error exception is thrown). Each
     * element of the parents vector gives the source of the value to
     * be copied.  A node may appear several times in this vector.
     *
     * @param offsets Vector of offsets. The offsets vector must have
     * the same size as the value array of the AggNode to be
     * constructed (or a length_error exception is thrown).  Each
     * element gives the element of the value vector of the
     * corresponding parent node to be copied. If the offset is greater
     * than the length of the corresponding parent, an out_of_range
     * exception is thrown.
     *
     * @exception length_error out_of_range
     */
    AggNode(std::vector<unsigned long> const &dim,
	    unsigned int nchain,
	    std::vector<Node const *> const &parents, 
	    std::vector<unsigned long> const &offsets);
    ~AggNode();
    /**
     * Copies values from parents.
     */
    void deterministicSample(unsigned int chain);
    /**
     * An aggregate node is discrete valued if all of its parents are.
     */
    bool isDiscreteValued() const;
    /**
     * Aggregate nodes are closed under all classes and are always fixed.
     */
    bool isClosed(std::set<Node const *> const &ancestors, 
		  ClosedFuncClass fc, bool fixed) const;
    /**
     * An AggNode places no restrictions on its parents' values. Therefore
     * this function always returns true.
     */
    bool checkParentValues(unsigned int chain) const;
    /**
     * An aggregate node is named after its first and last parents
     */
    std::string deparse(std::vector<std::string> const &parents) const;
    /**
     * Returns the vector of offsets
     */
    std::vector<unsigned long> const &offsets() const;
    //DeterministicNode *clone(std::vector<Node const *> const &parents) const;
};

} /* namespace jags */

#endif /* AGGREGATE_NODE_H */
