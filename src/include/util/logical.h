#ifndef LOGICAL_H_
#define LOGICAL_H_

#include <vector>
#include <algorithm>

namespace jags {

    /**
     * Tests whether all elements of the boolean vector "mask" are true.
     */
    inline bool allTrue (std::vector<bool> const &mask)
    {
	return std::find(mask.begin(), mask.end(), false) == mask.end();
    }
    
    /**
     * Tests whether any elements of the boolean vector "mask" are true.
     */
    inline bool anyTrue (std::vector<bool> const &mask)
    {
	return std::find(mask.begin(), mask.end(), true) != mask.end();
    }

    inline bool allFalse (std::vector<bool> const &mask)
    {
	return std::find(mask.begin(), mask.end(), true) == mask.end();
    }

    /**
     * Tests whether any elements of the boolean vector "mask" are true.
     */
    inline bool anyFalse (std::vector<bool> const &mask)
    {
	return std::find(mask.begin(), mask.end(), false) != mask.end();
    }

    /**
     * @short Gets a constant pointer to a unique boolean vector 
     *
     * This function creates a unique pointer to the requested vector,
     * avoiding redundant copies of the vector taking up memory.
     */
    std::vector<bool> const * getUnique(std::vector<bool> const &mask);

}

#endif /* LOGICAL_H_ */
