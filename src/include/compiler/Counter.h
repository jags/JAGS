#ifndef COUNTER_H_
#define COUNTER_H_

#include <sarray/RangeIterator.h>

/**
 * @short indices for for loops
 *
 * The Counter class is used to represent counters in "for" loops in the
 * BUGS language. A Counter is a scalar RangeIterator.
 *
 * <pre>
 * for (i in 1:N) {
 * }
 * </pre>
 */
class Counter : public RangeIterator
{
public:
  /**
   * Create a counter
   * @param range Scalar range. i.e., Range having ndim(false) == 1
   */
  Counter(Range const &range);
  /** Increments the value of the counter by 1 */
  Counter &next();
};

#endif /* COUNTER_H_ */

