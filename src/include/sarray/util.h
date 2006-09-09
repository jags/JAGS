#ifndef UTIL_H_
#define UTIL_H_

/**
 * @short Some low-level utility functions.
 */
#include <sarray/Index.h>
#include <vector>
#include <algorithm>

/**
 * Tests whether all elements of the boolean vector "mask" are true.
 */
inline bool allTrue (std::vector<bool> const &mask)
{
   return find(mask.begin(), mask.end(), false) == mask.end();
}

/**
 * Tests whether any elements of the boolean vector "mask" are true.
 */
inline bool anyTrue (std::vector<bool> const &mask)
{
   return find(mask.begin(), mask.end(), true) != mask.end();
}

/**
 * Tests whether the dimension represented by the Index "dim"
 * corresponds to a scalar quantity.
 *
 * Note that isScalar, isVector, and isArray are mutually exclusive.
 */
inline bool isScalar(Index const &dim)
{
   return dim.size() == 1 && dim[0] == 1;
}

/**
 * Tests whether the dimension represented by the Index "dim"
 * corresponds to a vector quantity.

 * Note that isScalar, isVector, and isArray are mutually exclusive.
 */
inline bool isVector(Index const &dim)
{
   return dim.size() == 1 && dim[0] > 1;
}

/**
 * Tests whether the dimension represented by the Index "dim"
 * corresponds to an array.
 *
 * Note that isScalar, isVector, and isArray are mutually exclusive.
 */
inline bool isArray(Index const &dim)
{
   return dim.size() > 1;
}

/**
 * Tests whether the dimension represented by the Index "dim"
 * corresponds to a two-dimensional array.
 */
inline bool isMatrix(Index const &dim)
{
   return dim.size() == 2;
}

/**
 * Tests whether the dimension represented by the Index "dim"
 * corresponds to a matrix with the same number of rows as columns.
 */
inline bool isSquareMatrix(Index const &dim)
{
   return dim.size() == 2 && dim[0] == dim[1];
}

#endif /* UTIL_H_ */
