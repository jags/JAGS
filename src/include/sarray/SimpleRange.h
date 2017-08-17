#ifndef SIMPLE_RANGE_H_
#define SIMPLE_RANGE_H_

#include <sarray/Range.h>

#include <string>
#include <vector>

namespace jags {

    /**
     * @short Represents a contiguous ordered collection of array indices 
     *
     * SimpleRange is used for a Range consisting of a contiguous
     * block of indices. It is uniquely defined by its lower and upper
     * bounds and contains every index in between.
     */
    class SimpleRange : public Range {
      public:
	/**
	 * Default constructor which constructs a NULL range, with
	 * zero-length upper and lower limits.
	 */
	SimpleRange();
	/**
	 * Destructor
	 */
	~SimpleRange();
	/**
	 * Constructs a simple range given lower and upper limits A
	 * logic_error is thrown if these are of different lengths
	 *
	 * @param lower Lower limits. 
	 *
	 * @param upper Upper limits. A range_error is thrown if any
	 *              element of upper is smaller than the corresponding
	 *              element of lower.
	 *
	 * @exception range_error
	 */
	SimpleRange(std::vector<int> const &lower, 
		    std::vector<int> const &upper);
	/**
	 * Constructs a scalar range from an index. The upper and lower
	 * bounds are both equal to the supplied index.  
	 */
	SimpleRange(std::vector<int> const &index);
	/**
	 * Constructs a range from a dimension. For each index, the lower
	 * limit is 1 and  the upper limit is the corresponding element of
	 * dim (cast to a signed int).
	 *
	 * This constructor should not be confused with the constructor
	 * that creates a scalar range from a vector of signed integers.
	 */
	SimpleRange(std::vector<unsigned int> const &dim);
	/**
	 * Equality operator. A SimpleRange is uniquely defined by
	 * its lower and upper limits, so testing for equality is
	 * much quicker than for a Range object.
	 */
	bool operator==(SimpleRange const &range) const;
	/**
	 * Inequality operator 
	 */
	bool operator!=(SimpleRange const &range) const;
	/*
	 * Indicates whether the supplied index is contained
	 * inside this SimpleRange.
	 bool contains(std::vector<int> const &test_index) const;
	*/
	/**
	 * Indicates whether the test SimpleRange is completely
	 * contained inside this SimpleRange.
	 *
	 * @param test_range Test SimpleRange
	 */
	bool contains(SimpleRange const &test_range) const;
	/**
	 * Indicates whether the test Range is completely contained
	 * inside this SimpleRange.
	 *
	 * @param test_range Test Range
	 */
	bool contains(Range const &test_range) const;
	/**
	 * The inverse of leftIndex. The left offset is the position
	 * of the given index in the sequence when all indices in the
	 * Range are written out in column major order (i.e. with the
	 * left index moving fastest).
	 *
	 * @param index Index vector to convert to offset. An out_of_range
	 * exception is thrown if the index is not contained in the range.
	 *
	 * @exception out_of_range
	 */
	unsigned int leftOffset(std::vector<int> const &index) const;
	/**
	 * The inverse of rightIndex. The right offset is the position
	 * of the given index in the sequence when all indices in the
	 * range are given in row major order (i.e. with the right
	 * index moving fastest).
	 *
	 * @param index Index vector to convert to offset. An out_of_range
	 * exception is thrown if the index is not contained in the range.  
	 *
	 * @exception out_of_range
	 */
	unsigned int rightOffset(std::vector<int> const &index) const;
	/**
	 * Less than operator. It defines the same ordering as
	 * Range#operator< among SimpleRange objects but the
	 * calculations are more efficient as a SimpleRange is
	 * uniquely defined by its first and last elements.
	 */
	bool operator<(SimpleRange const &rhs) const;
	/**
	 * The lower bound of the Range (an alias for Range#first)
	 */
	inline std::vector<int> const & lower() const { return first(); }
	/**
	 * The upper bound of the Range (an alias for Range#last)
	 */
	inline std::vector<int> const & upper() const { return last(); }
    };

    /**
     * Returns a string containing a BUGS language representation of
     * the given range: e.g. a simple range with lower limit (1,2,3)
     * and upper limit (3,3,3) will be printed as "[1:3,2:3,3]"
     */
    std::string print(SimpleRange const &range);
    
} /* namespace jags */

#endif /* SIMPLE_RANGE_H_ */
