#include <config.h>
#include <sarray/SArray.h>
#include "Logit.h"

#include <cmath>

using std::vector;
using std::log;

Logit::Logit ():ScalarFunc ("logit", 1)
{
}

double Logit::eval (vector <SArray const *> const &args) const
{
  double arg = *args[0]->value ();
  return log (arg / (1 - arg));
}

bool Logit::checkParameterValue (vector <SArray const *> const &args) const
{
  double arg = *args[0]->value ();
  return (arg >= 0 && arg <= 1);
}
