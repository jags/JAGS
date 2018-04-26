#ifndef NODE_ARRAY_SUBSET_H_
#define NODE_ARRAY_SUBSET_H_

#include <sarray/Range.h>

#include <vector>

namespace jags {

    class Node;
    class NodeArray;
    class Range;
    
    /**
     * @short Subset of a NodeArray
     * 
     */
    class NodeArraySubset {
	std::vector<unsigned long> _dim;
	unsigned int _nchain;
	std::vector<Node *> _node_pointers;
	std::vector<unsigned long> _offsets;
      public:
	/**
	 * Constructor. Creates a NodeArraySubset from a NodeArray
	 * and a given range
	 */
	NodeArraySubset(NodeArray const *array, Range const &range);
	/**
	 * Returns the values of the nodes in the range covered by the
	 * NodeArraySubset in column major order.
	 *
	 * @param chain Index number of chain to read.
	 */
	std::vector<double> value(unsigned int chain) const;
	/**
	 * Returns the dimension of the subset
	 */
	std::vector<unsigned long> const &dim() const;
	/**
	 * Returns a vector containing the nodes that contribute
	 * values to the NodeArraySubset. Repeated values are removed.
	 */
	std::vector<Node const *> nodes() const;
	/**
	 * Returns a vector containing all the nodes that contribute
	 * values to the NodeArraySubset, including repeated values.
	 */
	std::vector<Node const *> allnodes() const;
	/**
	 * Returns the number of chains
	 */
	unsigned int nchain() const;
	/**
	 * Returns the number of values in the subset
	 */
	unsigned long length() const;
    };
    
} /* namespace jags */

#endif /* NODE_ARRAY_SUBSET_H */
