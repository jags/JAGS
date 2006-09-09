#ifndef SARRAY_H_
#define SARRAY_H_

#include <sarray/Range.h>

/** @short missing value */
extern const double JAGS_NA;

/**
 * An SArray represents a multi-dimensional array of double precision floats
 * @short multi-dimensional array
*/
class SArray 
{
  Range const _range;
  double *_value;
  bool  _fixed;
  bool _discrete;
  SArray &operator=(SArray const &rhs);
public:
  /**
   * Constructor for SArrays. 
   * 
   * On construction, the elements of the value array are all equal to
   * JAGS_NA.  
   * 
   * @param dim Dimension of SArray to be constructed
   */
  SArray(Index const &dim);
  /**
   * Copy constructor. 
   *
   * Note that the assignment operator of the SArray class is declared
   * private, and so cannot be used.  An SArray can only be copied by a
   * constructor. This ensures that, once constructed, an SArray cannot
   * change its length or dimension.
   *
   * A copy of a fixed SArray is not itself fixed, but is otherwise
   * identical to the original.
   */
  SArray(SArray const &orig);
  ~SArray();
  /**
   * Sets the value of an SArray. 
   *
   * It is an error to attempt to set the value of a fixed SArray.
   * A logic_error will be thrown in this case.
   *
   * @param value Array of values to be assigned.
   *
   * @param length Length of value argument. This must match the length
   * of the SArray or a length_error will be thrown.
   *
   * @exception logic_error length_error
   */
  void setValue(double const *value, unsigned long length);
  /**
   * Sets the value of a single element of SArray
   *
   * @param value Scalar value to be assigned
   *
   * @param offset Distance along the value array
   */
  void setValue(double value, unsigned long offset);
  /**
   * The value of the SArray in vector form. 
   *
   * Values are given in column-major order (i.e. with
   * the left hand index moving fastest), like the S language and Fortran. 
   *
   * @return A pointer to the beginning of the value array.
   */
  inline double const *value() const;
  /** 
   * Sets the SArray to be fixed (true) or not (false).
   *
   * Fixed SArrays cannot change their value. It is an error to fix
   * an SArray that contains missing values.
   */
  void setFixed(bool fix);
  /**
   * Indicates whether the SArray is fixed, and so cannot change
   * its value. 
   */
  bool isFixed() const;
  /**
   * Determines whether the array is integer-valued (true) or real
   * valued (false).
   *
   * Currently, the integer-value flag of an SArray is purely
   * semantic: the storage type of the values is always double, and
   * the values are never checked to see if they are consistent with
   * the integer-value flag. This semantic information may be used, however,
   * in other parts of the library to see whether the value of the
   * SArray would be a valid subscript, for example.
   */
  void setDiscreteValued(bool flag);
  /**
   * Indicates whether the SArray should contain integer-values or
   * real values.
   */
  bool isDiscreteValued() const;
  /**
   * Returns the range associated with the SArray
   */
  Range const &range() const;
  /**
   * It is convenient to inline these functions so that an SArray
   * can be thought of as having some of the dimension attributes
   * of its associated range.
   */
  unsigned long length() const { return range().length(); }
  unsigned int ndim(bool drop) const { return range().ndim(drop); }
  Index const &dim(bool drop) const { return range().dim(drop); }
};


inline double const *SArray::value() const
{
  return _value;
}

#endif /* SARRAY_H_ */
