#include <config.h>
#include <sarray/SArray.h>
#include "Min.h"

#include <algorithm>

using std::min_element;
using std::min;
using std::vector;

Min::Min ()
  : ScalarFunc ("min", 0)
{
}

double Min::eval (vector <SArray const *> const &args) const
{
  double ans;
  for (unsigned int i = 0; i < args.size(); ++i) {
    double const *argi = args[i]->value ();
    long len = args[i]->length ();
    double mini = *min_element(argi, argi + len);
    if (i == 0) {
      ans = mini;
    }
    else {
      ans = min(ans, mini);
    }
  }
  return ans;
}

// Overriding ScalarFunc::checkParameterDim
bool Min::checkParameterDim (vector <Index> const &dims) const
{
  return true;
}

bool Min::isDiscreteValued(vector<bool> const &mask) const
{
    return count(mask.begin(), mask.end(), false) == 0;
}
