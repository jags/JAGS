#ifndef FUNC_LOGIT_H_
#define FUNC_LOGIT_H_

#include <function/ScalarFunc.h>

/**
 * @short Logistic function
 * @see ILogit
 * <pre>
 * y <- logit(p)
 * y = log(p/(1-p)) for 0 < p < 1
 * </pre>
 */
class Logit:public ScalarFunc
{
public:
  Logit ();
  double eval (std::vector <SArray const *> const &args) const;
  bool checkParameterValue (std::vector <SArray const *> const &args) const;
};

#endif /* FUNC_LOGIT_H_ */
