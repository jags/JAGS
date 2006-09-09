#ifndef FUNC_PROBIT_H_
#define FUNC_PROBIT_H_

#include <function/ScalarFunc.h>

/**
 * @short Normal quantile function
 * @see Phi
 * <pre>
 * y <- probit(p)
 * Prob(Y < y) = p for Y ~ N(0,1); 0 < p < 1
 * </pre>
 */
class Probit:public ScalarFunc
{
public:
  Probit ();
  double eval (std::vector <SArray const *> const &args) const;
  bool checkParameterValue (std::vector <SArray const *> const &args) const;
};

#endif /* FUNC_PROBIT_H_ */
