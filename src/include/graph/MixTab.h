#ifndef MIX_MAP_H_
#define MIX_MAP_H_

#include <sarray/SimpleRange.h>

#include <vector>
#include <map>

namespace jags {

    class Node;

    /**
     * @short Helper class for MixtureNode
     *
     * Evaluation of a MixtureNode uses a MixTab, which takes the
     * current map value (a vector of integers) and returns a
     * pointer to the parent node from which the MixtureNode takes its
     * current value.
     */
    class MixTab {
	const SimpleRange _range;
	std::vector<Node const *> _nodes;
      public:
	/**
	 * Constructs a MixTab from a MixMap
	 */
	MixTab(std::map<std::vector<unsigned long>, Node const *> const &mixmap);
	/** 
	 * Returns a pointer to the node corresponding to the given
	 * index.  If there is no node matching the index, a NULL
	 * pointer is returned.
	 */
	Node const * getNode(std::vector<unsigned long> const &index) const;
	/**
	 * Returns the range covered by the indices
	 */
	Range const &range() const;
    };

} /* namespace jags */

#endif /* MIX_MAP_H_ */
