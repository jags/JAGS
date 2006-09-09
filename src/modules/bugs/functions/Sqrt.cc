#include <config.h>
#include <sarray/SArray.h>
#include "Sqrt.h"

#include <cmath>

using std::vector;
using std::sqrt;

Sqrt::Sqrt ():ScalarFunc ("sqrt", 1)
{
}

double
Sqrt::eval (vector <SArray const *> const &args) const
{
  return sqrt(args[0]->value()[0]);
}

bool
Sqrt::checkParameterValue (vector <SArray const *> const &args) const
{
  return *args[0]->value() >= 0;
}
