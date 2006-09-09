#ifndef SCALAR_FUNC_H_
#define SCALAR_FUNC_H_

#include <function/Function.h>

/** 
 * ScalarFunc provides a base class for scalar-valued functions
 * @short Scalar functions
 */
class ScalarFunc : public Function
{
public:
  ScalarFunc (std::string const &name, unsigned int npar);
  /** 
   * Implementation of the virtual function Function::evaluate, which
   * calls ScalarFunc##evalScalar.
   */
  void evaluate (SArray & value,
		 std::vector <SArray const *> const &args) const;
  /**
   * Scalar functions need to implement this simplified function,
   * which returns a double, instead of Function##evaluate.
   */
  virtual double eval (std::vector <SArray const *> const &args) const =
    0;
  /**
   * Returns Index(1) since the result of a scalar function should
   * always be scalar
   */
  Index dim (std::vector <Index> const &args) const;
  /**
   * Most functions with a scalar value also take scalar arguments.
   * For convenience, a default implementation of Function##checkParameterDim
   * is provided which checks that all arguments are scalar.  This 
   *  must be overridden if the function takes any vector arguments
   */
  bool checkParameterDim (std::vector <Index> const &args) const;
};

#endif /* SCALAR_FUNC_H_ */
