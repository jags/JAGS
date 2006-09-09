#include <config.h>
#include <sarray/SArray.h>
#include "Pow.h"

#include <cmath>
#include <cfloat>

using std::vector;
using std::fabs;
using std::pow;

Pow::Pow ()
  : ScalarFunc ("pow", 2)
{
}

double Pow::eval (vector <SArray const *> const &args) const
{
  double arg1 = *(args[0]->value ());
  double arg2 = *(args[1]->value ());
  return pow (arg1, arg2);
}

bool Pow::checkParameterValue (vector <SArray const *> const &args) const
{
  if (*(args[0]->value()) >= 0) {
    return true;
  }
  else {
    double arg2 = *(args[1]->value ());
    int iarg2 = static_cast<int>(arg2 + DBL_EPSILON);
    return fabs(arg2 - iarg2) < DBL_EPSILON;
  }
}
