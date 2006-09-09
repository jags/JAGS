#include <config.h>
#include <sarray/SArray.h>
#include <sarray/Range.h>
#include "InProd.h"

using std::vector;

InProd::InProd () : ScalarFunc ("inprod", 2)
{
}

double
InProd::eval (vector <SArray const *> const &args) const
{
  double const *arg1 = args[0]->value ();
  double const *arg2 = args[1]->value ();
  long len = args[0]->length ();
  double svalue = 0;
  for (long i = 0; i < len; i++)
    {
      svalue += arg1[i] * arg2[i];
    }
  return svalue;
}

bool InProd::checkParameterDim (vector <Index> const &dims) const
{
    return (dims[0] == dims[1]);
}

bool InProd::isDiscreteValued(std::vector<bool> const &mask) const
{
    return count(mask.begin(), mask.end(), false) == 0;
}

bool InProd::isLinear(vector<bool> const &mask, vector<bool> const &fix) const
{
  //Test for quadratic term
  if (mask[0] && mask[1])
     return false;

  if (!fix.empty()) {
     if ((!mask[0] && !fix[0]) || (!mask[1] && !fix[1]))
        return false;
  }

  return true;
}

bool InProd::isScale(unsigned int index, vector<bool> const &fix) const
{
  if (fix.empty()) {
     return true;
  }
  else {
     if ((index == 0 && !fix[1]) || (index == 1 && !fix[0])) {
        return false;
     }
  }

  return true;
}


