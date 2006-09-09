#include <config.h>
#include <sarray/SArray.h>
#include <sarray/util.h>
#include "Add.h"

using std::vector;

Add::Add ():ScalarFunc ("+", 0)
{
}

double
Add::eval (vector <SArray const *> const &args) const
{
  double out = *args[0]->value();
  unsigned int n = args.size();
  for (unsigned int i = 1; i < n; ++i) {
    out += *args[i]->value();
  }
  return out;
}

bool Add::isDiscreteValued(vector<bool> const &mask) const
{
  return allTrue(mask);
}

bool Add::isLinear(vector<bool> const &mask, vector<bool> const &fix) const
{
    return true;
}

