#ifndef FUNC_STEP_H_
#define FUNC_STEP_H_

#include <function/ScalarFunc.h>

/**
 *@short Step function
 * <pre>
 * y <- step(x)
 * y = 1 if x >= 0
 *     0 if x <  0
 * </pre>
 */
class Step:public ScalarFunc
{
public:
  Step ();
  double eval (std::vector <SArray const *> const &args) const;
  bool isDiscreteValued(std::vector<bool> const &mask) const;
};

#endif /* FUNC_STEP_H_ */
