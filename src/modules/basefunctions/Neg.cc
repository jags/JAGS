#include <config.h>
#include <sarray/SArray.h>
#include <sarray/util.h>
#include "Neg.h"

#include <stdexcept>
#include <cmath>

using std::vector;
using std::logic_error;

Neg::Neg ()
  : ScalarFunc ("NEG", 1)
{
}

double
Neg::eval (vector <SArray const*> const &args) const
{
  double arg = *(args[0]->value ());
  return -arg;
}

bool Neg::isDiscreteValued(vector<bool> const &mask) const
{
  return allTrue(mask);
}

bool Neg::isLinear(vector<bool> const &mask, vector<bool> const &fix) const
{
    return true;
}

bool Neg::isScale(unsigned int index, vector<bool> const &fix) const
{
    return true;
}
