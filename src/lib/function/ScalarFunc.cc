#include <config.h>
#include <function/ScalarFunc.h>
#include <sarray/util.h>
#include <sarray/SArray.h>

using std::vector;
using std::string;

ScalarFunc::ScalarFunc (string const &name, unsigned int npar)
  : Function (name, npar)
{
}

void
ScalarFunc::evaluate (SArray & value, vector <SArray const *> const &args)
  const
{
  value.setValue(eval(args), 0);
}

Index ScalarFunc::dim (vector <Index> const &args) const
{
  return Index(1);
}

bool ScalarFunc::checkParameterDim (vector <Index> const &dims) const
{
  return count_if(dims.begin(), dims.end(), isScalar) == dims.size();
}
