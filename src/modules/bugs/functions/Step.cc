#include <config.h>
#include <sarray/SArray.h>
#include "Step.h"

using std::vector;

Step::Step ()
  : ScalarFunc ("step", 1)
{
}

double
Step::eval (vector <SArray const *> const &args) const
{
  return args[0]->value ()[0] >= 0 ? 1 : 0;
}

bool Step::isDiscreteValued(vector<bool> const &mask) const
{
  return true;
}
