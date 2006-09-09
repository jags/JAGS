#ifndef FUNC_ILOGIT_H_
#define FUNC_ILOGIT_H_

#include <function/InverseLinkFunc.h>

/**
 * @short Inverse of the logistic link function
 * @see Logit
 * <pre>
 * logit(y) <- a + b*x
 * y <- ilogit(a + b*x)
 * y = 1/(1 + exp(-a - b*x))
 * </pre>
 */
class ILogit:public InverseLinkFunc
{
public:
  ILogit ();
  double eval (std::vector <SArray const *> const &args) const;
  double link(double mu) const;
  double gradLink(double mu) const;
};

#endif /* FUNC_ILOGIT_H_ */
