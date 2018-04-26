#ifndef NODE_FACTORY_H_
#define NODE_FACTORY_H_

#include <vector>
#include <cfloat>

namespace jags {

    class Node;

    /**
     * @short Fuzzy less than function for doubles
     *
     * Gives an ordering for doubles allowing for some numerical
     * tolerance.  Values that are "nearly equal" are considered to be
     * equal.
     */
    inline bool lt(double value1, double value2)
    {
	return value1 < value2 - 16 * DBL_EPSILON;
    }
    
    /**
     * @short Fuzzy less than function for arrays of doubles
     *
     * Gives an ordering for two arrays of doubles of the same length.
     * Values that are numerically "nearly equal" are considered to
     * be equal.
     */
    bool lt(double const *value1, double const *value2, unsigned long length);

    /**
     * @short Fuzzy less than function for Nodes
     *
     * Gives an ordering for Nodes that is arbitrary, but reproducible
     * within the same process. This is all that is required by the
     * various Node factories.
     * 
     * In this ordering fixed nodes come before non-fixed nodes. Fixed
     * nodes are ordered within themselves using the fuzzy less than
     * operator for their values. Non-fixed nodes are sorted by their
     * pointer values (which gives an arbitrary, but strict ordering).
     */
    bool lt(Node const *node1, Node const *node2);

    /**
     * @short Fuzzy less than function for vectors of Nodes
     *
     * Extends the fuzzy less than operator for Nodes to vectors of
     * Nodes. Vectors are first sorted by length, then an element-wise
     * comparison of each Node in the vector.
     */
    bool lt(std::vector<Node const *> const &par1, 
	    std::vector<Node const *> const &par2);
    

    /**
     * Template that constructs a function object from any of the
     * above less than functions.
     */
    template<typename T> struct fuzzy_less
    {
	bool operator()(T const &arg1, T const &arg2) const {
	    return lt(arg1, arg2);
	}
    };

} /* namespace jags */

#endif /* NODE_FACTORY_H_ */
