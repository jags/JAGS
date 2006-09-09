#ifndef FUNC_LOGGAM_H_
#define FUNC_LOGGAM_H_

#include <function/ScalarFunc.h>

/**
 * @short Log gamma function
 * <pre>
 * y <- loggam(x)
 * y = log(gamma(x)) for x > 0
 * </pre>
 */
class LogGam:public ScalarFunc
{
public:
  LogGam ();
  double eval (std::vector <SArray const *> const &args) const;
  bool checkParameterValue (std::vector <SArray const *> const &args) const;
};

#endif /* FUNC_LOGGAM_H_ */
