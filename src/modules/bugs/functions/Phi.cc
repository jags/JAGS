#include <config.h>
#include <sarray/SArray.h>
#include "Phi.h"

#include <Rmath.h>

using std::vector;

Phi::Phi ():InverseLinkFunc ("phi", "probit")
{
}

double Phi::eval (vector <SArray const *> const &args) const
{
  double q = *(args[0]->value ());
  if (!R_FINITE (q)) {
    return q > 0 ? 1 : 0;
  }
  double p = pnorm5 (q, 0, 1, 0, 0);
  if (p == 0)
    return DBL_EPSILON;
  else if (p == 1)
    return 1 - DBL_EPSILON;
  else
    return p;
}

double Phi::link(double mu) const
{
  return pnorm (mu, 0, 1, 0, 0);
}

double Phi::gradLink(double mu) const
{
  return dnorm (mu, 0, 1, 0);
}
