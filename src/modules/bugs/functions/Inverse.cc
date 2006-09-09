#include <config.h>
#include "matrix.h"
#include <sarray/SArray.h>
#include <sarray/util.h>
#include "Inverse.h"

using std::vector;

Inverse::Inverse ():Function ("inverse", 1)
{
}

void
Inverse::evaluate (SArray & value, vector <SArray const *> const &args) const
{
    int n = args[0]->dim(true)[0];
    int N = args[0]->length();
    double *y = new double[N];
    inverse (y, args[0]->value(), n, false);
    value.setValue (y, N);
    delete y;
}

bool Inverse::checkParameterDim (vector <Index> const &dims) const
{
    return isSquareMatrix(dims[0]);
}

Index Inverse::dim (vector <Index> const &dims) const
{
    return dims[0];
}
