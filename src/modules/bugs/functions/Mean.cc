#include <config.h>
#include <sarray/SArray.h>
#include "Mean.h"

using std::vector;

Mean::Mean ()
    : ScalarFunc ("mean", 1)
{
}

double
Mean::eval (vector <SArray const*> const &args) const
{
    double const *arg1 = args[0]->value ();
    long len = args[0]->length ();
    double svalue = 0;
    for (long i = 0; i < len; i++)
	{
	    svalue += arg1[i];
	}
    svalue /= len;
    return svalue;
}

bool Mean::checkParameterDim (vector <Index> const &dims) const
{
    return true;
}

bool Mean::isLinear(vector<bool> const &mask, vector<bool> const &fix) const
{
    return true;
}

bool Mean::isScale(unsigned int index, vector<bool> const &fix) const
{
    return true;
}

