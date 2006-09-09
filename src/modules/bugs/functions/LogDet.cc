#include <config.h>

#include "LogDet.h"
#include "matrix.h"

#include <sarray/util.h>
#include <sarray/SArray.h>

using std::vector;

LogDet::LogDet ()
  : ScalarFunc ("logdet", 1)
{
}

double LogDet::eval (vector <SArray const *> const &args) const
{
  return logdet (args[0]->value(), args[0]->dim(true)[0]);
}

bool LogDet::checkParameterDim (vector <Index> const &dims) const
{
  return isSquareMatrix(dims[0]);
}

/* FIXME: we need a checkParameterValue function */

