#include <config.h>
#include <sarray/util.h>
#include <sarray/SArray.h>
#include "Rank.h"

#include <algorithm>

using std::vector;
using std::stable_sort;

bool lt_doubleptr (double const *arg1, double const *arg2) {
  return *arg1 < *arg2;
}

Rank::Rank ()
  : Function ("rank", 1)
{
}

void Rank::evaluate (SArray & value, vector <SArray const *> const &args) const
{
  double const *arg = args[0]->value ();
  long len = args[0]->length ();

  //Create a vector of pointers to the elements of arg and sort it
  double const **argptrs = new double const *[len];
  for (long i = 0; i < len; ++i) {
    argptrs[i] = arg + i;
  }
  stable_sort(argptrs, argptrs + len, lt_doubleptr);

  //Ranks can be inferred from the sorted vector of pointers
  double *ranks = new double[len];
  for (long i = 0; i < len; ++i) {
    ranks[argptrs[i] - arg] = i + 1;
  }
  delete [] argptrs;
  
  value.setValue(ranks, len);
  delete [] ranks;
}

Index Rank::dim (vector <Index> const &dims) const
{
  return dims[0];
}

bool Rank::checkParameterDim (vector<Index> const &dims) const
{
    return isVector(dims[0]) || isScalar(dims[0]);
}

bool Rank::isDiscreteValued(vector<bool> const &mask) const
{
  return true;
}
