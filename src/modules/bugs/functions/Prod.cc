#include <config.h>
#include <sarray/util.h>
#include "Prod.h"

using std::vector;

Prod::Prod () : Function("prod", 1)
{
}

void Prod::evaluate(double *x, vector <double const *> const &args,
	      vector<vector<unsigned int> > const &dims) const
{
    double value = args[0][0];
    unsigned int len = product(dims[0]);
    for (unsigned int i = 1; i < len; ++i) {
	value *= args[0][i];
    }
    *x = value;
}

bool Prod::checkParameterDim (vector<vector<unsigned int> > const &args) const
{
    return true;
}

bool Prod::isDiscreteValued(vector<bool> const &mask) const
{
    return allTrue(mask);
}