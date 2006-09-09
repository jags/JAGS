#include <config.h>
#include <sarray/SArray.h>
#include "Sum.h"

#include <stdexcept>

using std::vector;
using std::logic_error;

Sum::Sum () : ScalarFunc ("sum", 1)
{
}

double
Sum::eval (vector <SArray const *> const &args) const
{
  double const *arg1 = args[0]->value ();
  long len = args[0]->length ();
  double value = 0;
  for (long i = 0; i < len; i++)
    {
      value += arg1[i];
    }
  return value;
}

bool Sum::checkParameterDim (vector <Index> const &args) const
{
  return true;
}

bool Sum::isDiscreteValued(vector<bool> const &mask) const
{
    return count(mask.begin(), mask.end(), false) == 0;
}

bool Sum::isLinear(vector<bool> const &mask, vector<bool> const &fix) const
{
  return true;
}

bool Sum::isScale(unsigned int index, vector<bool> const &fix) const
{
  return true;
}
