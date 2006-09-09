#ifndef RANGE_IMP_H_
#define RANGE_IMP_H_

#include <sarray/Index.h>

/**
 * A RangeImp implements the interface described by the Range class.
 * All of the member functions of Range are constant, so identical
 * Range objects can share memory.
 * 
 * The interface is documented in the Range class.
 * 
 * @see Range
 */
class RangeImp {
  Index _lower, _upper, _dim, _dim_dropped;
  unsigned long _length;
  /** 
   * Forbid copying. We only need one instance of a RangeImp object
   * with the same upper and lower limits 
   */
  RangeImp(RangeImp const &rhs);
  RangeImp &operator=(RangeImp const &rhs);
 public:
  /** @see Range##Range */
  RangeImp(Index const &lower, Index const &upper);
  /** @see Range##length */
  unsigned long length() const;
  /** @see Range##contains */
  bool contains(RangeImp const &range) const;
  /** @see Range##overlaps */
  bool overlaps(RangeImp const &range) const;
  /** @see Range##leftOffset */
  unsigned long leftOffset(Index const &index) const;
  /** @see Range##leftIndex */
  Index leftIndex(unsigned long offset) const;
  /** @see Range##rightOffset */
  unsigned long rightOffset(Index const &index) const;
  /** @see Range##rightIndex */
  Index rightIndex(unsigned long offset) const;
  /** @see Range##Index */
  Index const & dim(bool drop) const;
  /** @see Range##ndim */
  unsigned int ndim(bool drop) const;
  /** @see Range##lower */
  Index const & lower() const;
  /** @see Range##upper */
  Index const & upper() const;
  /** @see Range##operator< */
  bool operator<(RangeImp const &rhs) const;
};


#endif /* RANGE_IMP_H_ */
