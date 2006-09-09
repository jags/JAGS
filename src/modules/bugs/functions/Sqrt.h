#ifndef FUNC_SQRT_H_
#define FUNC_SQRT_H_

#include <function/ScalarFunc.h>

/**
 * @short Square root function
 * <pre>
 * sigma <- sqrt(tau)
 * tau = sigma * sigma; tau > 0
 * </pre>
 */
class Sqrt : public ScalarFunc
{
public:
  Sqrt ();
  double eval (std::vector <SArray const *> const &args) const;
  bool checkParameterValue (std::vector <SArray const *> const &args) const;
};

#endif /* FUNC_SQRT_H_ */
