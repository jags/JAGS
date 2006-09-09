#include <config.h>
#include <sarray/SArray.h>
#include <sarray/util.h>
#include "Subtract.h"

#include <stdexcept>

using std::vector;
using std::logic_error;

Subtract::Subtract ()
  : ScalarFunc ("-", 2)
{
}

double
Subtract::eval (vector <SArray const *> const &args) const
{
  double val1 = *(args[0]->value ());
  double val2 = *(args[1]->value ());
  return val1 - val2;
}

bool Subtract::isDiscreteValued(vector<bool> const &mask) const
{
  return allTrue(mask);
}

bool 
Subtract::isLinear(vector<bool> const &mask, vector<bool> const &fix) const
{
    return true;
}
