#ifndef FUNC_LOG_H_
#define FUNC_LOG_H_

#include <function/ScalarFunc.h>

/**
 * @short Log function
 * @see Exp
 * <pre>
 * y <- log(x)
 * y = log(x) for x > 0
 * </pre>
 */
class Log : public ScalarFunc
{
public:
  Log ();
  double eval (std::vector <SArray const *> const &args) const;
  bool checkParameterValue (std::vector <SArray const *> const &args) const;
};

#endif /* FUNC_LOG_H_ */
