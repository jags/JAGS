#include <config.h>
#include <sarray/SArray.h>
#include "Equals.h"

#include <cmath>
#include <cfloat>

using std::vector;
using std::fabs;

Equals::Equals ():ScalarFunc ("equals", 2)
{
}

double Equals::eval (vector <SArray const *> const &args) const
{
  return (fabs (*args[0]->value() - *args[1]->value()) < DBL_EPSILON);
}

bool Equals::isDiscreteValued(std::vector<bool> const &mask) const
{
  return true;
}
