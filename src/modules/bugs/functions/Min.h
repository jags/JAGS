#ifndef FUNC_MIN_H_
#define FUNC_MIN_H_

#include <function/ScalarFunc.h>

/**
 * @short Minimum element of a vector
 * @see Max
 */
class Min: public ScalarFunc
{
public:
  Min ();
  double eval (std::vector <SArray const *> const &args) const;
  bool checkParameterDim (std::vector <Index> const &dims) const;
  Index dim (std::vector <SArray const *> const &args) const;
  /** Returns true if all arguments are true */
  bool isDiscreteValued(std::vector<bool> const &mask) const;  
};

#endif /* FUNC_MIN_H_ */
