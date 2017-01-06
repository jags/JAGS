#ifndef RANGE_H_
#define RANGE_H_

#include <string>
#include <vector>

namespace jags {

/**
 * @short Represents a collection of array indices 
 *
 * A Range object represents a subset of indices in a multi
 * dimensional array. For example, A[1:2,4,4:6] is a subset of the
 * 3-dimensional array A defined by the range [1:2,4,4:6]. This range
 * is itself defined as the outer product of 3 integer vectors:
 * c(1,2), 4, and c(4,5,6) in the first, second, and third dimensions,
 * respectively.
 *
 * A Range may contain indices in any order and may include repeat
 * indices. A Range is "simple" if its component vectors consist of
 * contiguous elements in increasing order. For example, the Range
 * given above is Simple.  However, the Range [c(2,1,1,3,4)] is not
 * simple. The SimpleRange class is used for a more efficient
 * calculations on simple Ranges.
 *
 * @see SimpleRange
 * 
 */
class Range {
    std::vector<std::vector<int> > _scope;
  protected:
    std::vector<unsigned int> _dim, _dim_dropped;
    std::vector<int> _first, _last;
    unsigned int _length;
  public:
    /**
     * Default constructor which creates a NULL range of zero length.
     */
    Range();
    /**
     * Constructs a Range.
     *
     * @param scope Vector of index vectors. The scope parameter is of
     * length equal to the number of dimensions. Element i of the
     * scope vector is the vector of indices for dimension i. Each
     * element of scope must be non-empty (i.e. must contain at least
     * one index) or a logic_error is thrown.
     *
     * @exception logic_error
     */
    Range(std::vector<std::vector<int> > const &scope);
    /**
     * Virtual destructor
     */
    virtual ~Range();
    /**
     * Equality operator
     */
    bool operator==(Range const &range) const;
    /**
     * Inequality operator 
     */
    bool operator!=(Range const &range) const;
    /**
     * Length of the range. This is the number of indices that are
     * contained in the range. For example, the range [1:2,4,4:5]
     * contains 4 elements.
     */
    unsigned int length() const;
    /**
     * Returns the index in position n when the indices in the Range
     * are put in column-major order (i.e. with the left hand index
     * moving fastest).
     *
     * @see RangeIterator
     */
    std::vector<int> leftIndex(unsigned int n) const;
    /**
     * Returns the index in position n when the indices in the Raneg
     * are put in row-major order (i.e. with the right hand index
     * moving fastest).
     *
     * @see RangeIterator
     */
    std::vector<int> rightIndex(unsigned int n) const;
    /**
     * Dimension of the range. The range [1:4,2,3:5] has dimension
     * (4,1,3) if drop==false and (4,3) if drop==true. Dropping of
     * dimensions follows the S language convention.
     *
     * @param drop Should dimensions of size 1 be dropped? 
     */
    std::vector<unsigned int> const &dim(bool drop) const;
    /**
     * Number of dimensions covered by the Range. The range [1:4, 2,
     * 3:5] has 3 dimensions if drop==false and 2 dimensions if
     * drop==true.
     *
     * @param drop Should dimensions of size 1 be counted?
     */
    unsigned int ndim(bool drop) const;
    /**
     * The first element of the Range (in either row-major or
     * column-major order).
     */
    std::vector<int> const & first() const;
    /**
     * The last element of the Range (in either row-major or 
     * column-major order).
     */
    std::vector<int> const & last() const;
    /**
     * Less than operator that gives a unique ordering of ranges
     * Ranges are first sorted by the first element, then the last
     * element and, if they are still equivalent by these criteria, a
     * lexicographic ordering of the scope.
     */
    bool operator<(Range const &rhs) const;
    /**
     * Returns the scope vector used to construct the Range. 
     */
    std::vector<std::vector<int> > const &scope() const;
};

/**
 * Tests for NULL ranges.
 */
inline bool isNULL(Range const &range) { return range.length() == 0; }

    /**
     * Returns a string containing a BUGS language representation of
     * the given range.
     *
     * Simple ranges, consisting of contiguous elements in increasing
     * order, are represented as in R: e.g. a range with lower limit
     * (1,2,3) and upper limit (3,3,3) will be printed as
     * "[1:3,2:3,3].
     *
     * Complex ranges, consisting of elements that are not contiguous
     * or not in order, are not represented exactly when printed. The
     * first and last indices in the range are given separated by
     * ellipses, e.g. the range c(3,7,4,2,2,5) is represented as
     * "[3...5]"
     */
    std::string print(Range const &range);

} /* namespace jags */

#endif /* RANGE_H_ */
