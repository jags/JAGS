#include <config.h>
#include <sarray/SArray.h>
#include "Probit.h"

#include <Rmath.h>

using std::vector;

Probit::Probit ()
  : ScalarFunc ("probit", 1)
{
}

double Probit::eval (vector <SArray const *> const &args) const
{
  double p = *(args[0]->value ());
  return qnorm5 (p, 0, 1, 0, 0);
}

bool Probit::checkParameterValue (vector <SArray const *> const &args) const
{
  double p = *(args[0]->value ());
  return (p > 0 && p < 1);
}
