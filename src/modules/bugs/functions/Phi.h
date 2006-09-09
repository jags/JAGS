#ifndef FUNC_PHI_H_
#define FUNC_PHI_H_

#include <function/InverseLinkFunc.h>

/**
 * @short Standard normal distribution function
 * @see Probit
 * <pre>
 * p <- phi(x)
 * probit(p) <- x
 * p = Prob(X < x) where X ~ N(0,1)
 * </pre>
 */
class Phi:public InverseLinkFunc
{
public:
  Phi ();
  double eval (std::vector < SArray const *>const &args) const;
  double link(double mu) const;
  double gradLink(double mu) const;
};

#endif /* FUNC_PHI_H_ */
