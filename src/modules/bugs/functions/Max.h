#ifndef FUNC_MAX_H_
#define FUNC_MAX_H_

#include <function/ScalarFunc.h>

/**
 * @short Maximum element of a vector
 * @see Min
 */
class Max: public ScalarFunc
{
public:
  Max ();
  double eval (std::vector <SArray const *> const &args) const;
  bool checkParameterDim (std::vector <Index> const &dims) const;
  /** Returns true if all arguments are true */
  bool isDiscreteValued(std::vector<bool> const &mask) const;
};

#endif /* FUNC_MAX_H_ */
