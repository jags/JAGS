#ifndef FUNC_EQUALS_H_
#define FUNC_EQUALS_H_

#include <function/ScalarFunc.h>

/**
 * @short Test for equality
 * <pre>
 * y <- equals(a,b)
 * y = 1 if a==b
 *     0 if a!=b
 * </pre>
 */
class Equals:public ScalarFunc
{
public:
  Equals ();
  double eval (std::vector <SArray const *> const &args) const;
  /** Returns true */
  bool isDiscreteValued(std::vector<bool> const &mask) const;
};

#endif /* FUNC_EQUALS_H_ */
