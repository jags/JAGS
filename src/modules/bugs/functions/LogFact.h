#ifndef FUNC_LOGFACT_H_
#define FUNC_LOGFACT_H_

#include <function/ScalarFunc.h>

/**
 * @short Log factorial function
 * @see LogGamma
 * <pre>
 * y <- logfact(x)
 * y = x! 
 *   = log(gamma(x + 1)) for x > -1
 * </pre>
 */
class LogFact : public ScalarFunc
{
public:
  LogFact ();
  double eval (std::vector <SArray const *> const &args) const;
  bool checkParameterValue (std::vector <SArray const *> const &args) const;
};

#endif /* FUNC_LOGFACT_H_ */
