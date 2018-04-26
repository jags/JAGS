#ifndef RANGE_ITERATOR_H_
#define RANGE_ITERATOR_H_

#include <vector>

#include <sarray/Range.h>

namespace jags {

    /**
     * @short Mutable index that traverses a Range
     *
     * A RangeIterator is a numeric vector that traverses a given
     * Range in row- or column-major order. 
     *
     * @see Range
     */
    class RangeIterator : public std::vector<unsigned long> {
	std::vector<std::vector<unsigned long> > _scope;
	std::vector<unsigned long> _dim;
	std::vector<unsigned long> _index;
	unsigned int  _atend;
	//Forbid assignment
	RangeIterator &operator=(std::vector<unsigned long> const &);
      public:
	/**
	 * Constructor. The initial value of a RangeIterator is the
	 * beginning of the range.
	 *
	 * It is not possible to construct a RangeIterator if any of
	 * the dimensions of the Range are zero because there can be
	 * no index value corresponding to a zero dimension. In this
	 * case a logic_error is thrown.
	 *
	 * @param range Range to traverse
	 */
	RangeIterator(Range const &range);
	/**
	 * Goes to the next index in column-major order, (i.e. moving
	 * the left hand index fastest). If the RangeIterator reaches
	 * the end of the Range, then a call to nextLeft will move it
	 * to the beginning again.
	 * 
	 * @return reference to self after incrementation
	 * @see nextRight
	 */
	RangeIterator &nextLeft();
	/**
	 * Goes to the next index in row-major order (i.e. moving the
	 * right hand index fastest) but otherwise behaves like
	 * nextLeft.
	 *
	 * @return reference to self after incrementation
	 * @see nextLeft
	 */
	RangeIterator &nextRight();
	/**
	 * Returns a numeric counter of the number of times the
	 * RangeIterator has gone past the end of the Range (and
	 * returned to the beginning) in a call to nexLeft or
	 * nextRight. The initial value is zero.
	 */
	unsigned int atEnd() const;
    };

} /* namespace jags */

#endif /* RANGE_ITERATOR_H_ */

