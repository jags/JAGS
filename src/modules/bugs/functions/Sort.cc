#include <config.h>
#include <sarray/SArray.h>
#include "Sort.h"

#include <algorithm>
#include <stdexcept>

#include <sarray/util.h>

using std::vector;
using std::sort;
using std::logic_error;

Sort::Sort ()
  : Function ("sort", 1)
{
}

void Sort::evaluate (SArray & value, vector <SArray const *> const &args) const
{
  double const *arg1 = args[0]->value ();
  long len = args[0]->length ();

  double *copy = new double[len];
  for (long i = 0; i < len; ++i) {
    copy[i] = arg1[i];
  }
  sort(copy, copy + len);

  value.setValue(copy, len);
  delete [] copy;
}

Index Sort::dim (vector <Index> const &dims) const
{
    return dims[0];
}

bool Sort::checkParameterDim (vector<Index> const &dims) const
{
    return isVector(dims[0]) || isScalar(dims[0]);
}

bool Sort::isDiscreteValued(vector<bool> const &mask) const
{
  return count(mask.begin(), mask.end(), false) == 0;
}
