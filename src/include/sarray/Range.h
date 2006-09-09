#ifndef RANGE_H_
#define RANGE_H_

#include <sarray/Index.h>
#include <string>

class RangeImp;

/**
 * A Range object represents a range of indices used to take a subset
 * of a multi-dimensional array. For example, A[1:2,4,4:5] is a subset
 * of the array A defined by a range with lower boundary (1,4,4) and
 * upper boundary (2,4,5).
 *
 * All of the member functions of the Range class are constant, with
 * the exception of the assignment operator, so identical Range
 * objects (i.e. those with the same lower and upper bounds) can share
 * their memory.  In fact, the Range class is just a wrapper around
 * the RangeImp class, which actually implements the member functions.
 *
 * @short Range of indices 
 */
class Range {
  RangeImp const *_p;
 public:
  /**
   * Default constructor which constructs a NULL range, with zero-length
   * upper and lower limits.
   */
  Range();
  /**
   * Constructs a range based on given lower and upper limits
   * A logic_error is thrown if these are of different lengths
   * @param lower lower limits. A range_error is thrown if any
   *              of these are < 1
   * @param upper upper limits. A range_error is thrown if any
   *              element of upper is smaller than the corresponding
   *              element of lower.
   * @exception logic_error range_error
   */
  Range(Index const &lower, Index const &upper);
  /**
   * Constructs a scalar range from an index. The upper and lower
   * bounds are both equal to the supplied index.  This constructor
   * is used for implicit conversion of Index to Range objects, so
   * that an Index can be supplied to any function in place of a 
   * Range parameter.
   */
  Range(Index const &index);
  /** Copy constructor */
  Range (Range const &rhs);
  /** Destructor */
  ~Range();
  /** Assignment operator */
  Range &operator=(Range const &rhs);
  /** Equality operator */
  bool operator==(Range const &range) const;
  /** Inequality operator */
  bool operator!=(Range const &range) const;
  /**
   * Length of the range. This is the number of indices that are
   * contained in the range. For example, the range [1:2,4,4:5]
   * contains 4 elements.
   */
  unsigned long length() const;
  /**
   * Indicates whether the test range is completely contained inside this
   * range.
   *
   * @param test_range Test range, which must have the correct number of
   * dimensions, or an invalid_argument exception is thrown.  
   *
   * @exception invalid_argument
   */
  bool contains(Range const &test_range) const;
  /**
   * Returns the value of a RangeIterator after n iterations of nextLeft()
   * @see RangeIterator
   */
  Index leftIndex(long n) const;
  /**
   * The inverse of leftIndex. Returns the number of iterations of 
   * RangeIterator##nextLeft required to reach index.
   *
   * @param index Index to convert to offset. An out_of_range is thrown
   *              if the index is not contained in the range.
   * @see RangeIterator
   */
  long leftOffset(Index const &index) const;
  /**
   * Returns the value of a RangeIterator after n iterations of nextRight()
   * @see RangeIterator
   */
  Index rightIndex(long n) const;
  /**
   * The inverse of rightIndex. Returns the number of iterations of 
   * RangeIterator##nextRight required to reach index.
   *
   * @param index Index to convert to offset. An out_of_range is thrown
   *              if the index is not contained in the range.
   * @see RangeIterator
   */
  long rightOffset(Index const &index) const;
  /**
   * Dimension of the range. The range [1:4,2,3:5] has dimension
   * (4,1,3) if drop==false and (4,3) if drop==true. 
   *
   * @param drop Should dimensions of size 1 be dropped?
   */
  Index const &dim(bool drop) const;
  /**
   * Number of dimensions covered by the Range
   * @param drop Should dimensions of size 1 be counted?
   */
  unsigned int ndim(bool drop) const;
  /**
   * lower limit of range
   */
  Index const & lower() const;
  /**
   * upper limit of range
   */
  Index const & upper() const;
  /**
   * Lexicographic ordering based on lower and then upper limit
   * @see Index##operator<
   */
  bool operator<(Range const &rhs) const;
};

/**
 * Tests for NULL ranges.
 */
inline bool isNULL(Range const &range) { return range.length() == 0; }

/**
 * Returns a string containing a BUGS language representation of the
 * given range: e.g. a range with lower limit (1,2,3) and upper limit
 * (3,3,3) will be printed as "[1:3,2:3,3]"
 */
std::string print(Range const &range);

#endif /* RANGE_H_ */
