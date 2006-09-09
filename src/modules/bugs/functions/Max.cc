#include <config.h>
#include <sarray/SArray.h>
#include "Max.h"

#include <algorithm>

using std::vector;
using std::max_element;
using std::max;

Max::Max ()
  : ScalarFunc("max", 0)
{
}

double Max::eval (vector < SArray const *>const &args) const
{
  double ans;
  for (unsigned int i = 0; i < args.size(); ++i) {
    double const *argi = args[i]->value ();
    long len = args[i]->length ();
    double maxi = *max_element(argi, argi + len);
    if (i == 0) {
      ans = maxi;
    }
    else {
      ans = max(ans, maxi);
    }
  }
  return ans;
}

bool Max::checkParameterDim (vector <Index> const &dims) const
{
  return true;
}

bool Max::isDiscreteValued(vector<bool> const &mask) const
{
    return count(mask.begin(), mask.end(), false) == 0;
}
