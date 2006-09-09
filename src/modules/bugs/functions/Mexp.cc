#include <config.h>

#include "Mexp.h"
#include "matexp.h"

#include <sarray/util.h>
#include <sarray/SArray.h>

using std::vector;

Mexp::Mexp()
  : Function("mexp",1)
{
}

void
Mexp::evaluate (SArray &value, vector <SArray const *> const &args) const
{
  int n = args[0]->dim(true)[0];
  int N = args[0]->length();
  double *expmat = new double[N];
  MatrixExpPade(expmat, args[0]->value(), n, 1);
  value.setValue(expmat, N);
  delete [] expmat;
}

Index Mexp::dim (vector <Index> const &dims) const
{
    return dims[0];
}

bool Mexp::checkParameterDim (vector <Index> const &dims) const
{
    return isSquareMatrix(dims[0]);
}
