#include <config.h>
#include <util/logical.h>
#include <util/nainf.h>
#include <module/ModuleError.h>
#include "DSum.h"

#include <cfloat>
#include <cmath>

using std::vector;
using std::fabs;
using std::sqrt;

namespace bugs {

DSum::DSum()
    : ArrayDist("dsum", 0)
{
}

bool DSum::isDiscreteValued(vector<bool> const &mask) const
{
    return allTrue(mask);
}

double DSum::logDensity(double const *x, unsigned int length, PDFType type,
			vector<double const *> const &par,
			vector<vector<unsigned int> > const &dims,
			double const *lower, double const *upper) const
{
    const double tol = sqrt(DBL_EPSILON);
    for (unsigned int i = 0; i < length; ++i) {
	double s = x[i];
	for (unsigned int j = 0; j < par.size(); ++j) {
	    s -= par[j][i];
	}
	if (fabs(s) > tol) {
	    return JAGS_NEGINF;
	}
    }
    return 0;
}

void DSum::randomSample(double *x, unsigned int length,
			  vector<double const *> const &par, 
			  vector<vector<unsigned int> > const &dims,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
    for (unsigned int i = 0; i < length; ++i) {
	x[i] = 0;
	for (unsigned int j = 0; j < par.size(); ++j) {
	    x[i] += par[j][i];
	}
    }
}

void DSum::typicalValue(double *x, unsigned int length,
			  vector<double const *> const &par,
			  vector<vector<unsigned int> > const &dims,
			  double const *lower, double const *upper) const
{
    for (unsigned int i = 0; i < length; ++i) {
	x[i] = 0;
	for (unsigned int j = 0; j < par.size(); ++j) {
	    x[i] += par[j][i];
	}
    }
}

bool DSum::isSupportFixed(vector<bool> const &fixmask) const
{
    return allTrue(fixmask);
}

unsigned int DSum::df(vector<vector<unsigned int> > const &dims) const
{
    return 0;
}

bool DSum::checkParameterValue(vector<double const *> const &params,
			       vector<vector<unsigned int> > const &dims) const
{
    return true;
}

bool DSum::checkParameterDim (vector<vector<unsigned int> > const &dims) const
{
    for (unsigned int i = 1; i < dims.size(); ++i) {
	if (dims[i] != dims[0])
	    return false;
    }
    return true;
}

bool DSum::checkParameterDiscrete(vector<bool> const &mask) const
{
    for (unsigned int i = 1; i < mask.size(); ++i) {
	if (mask[i] != mask[0])
	    return false;
    }
    return true;
}

void DSum::support(double *lower, double *upper, unsigned int length,
		   vector<double const *> const &par,
		   vector<vector<unsigned int> > const &dims) const
{
    for (unsigned int i = 0; i < length; ++i) {
	lower[i] = 0;
	for (unsigned int j = 0; j < par.size(); ++j) {
	    lower[i] += par[j][i];
	}
	upper[i] = lower[i];
    }
}

vector<unsigned int> DSum::dim(vector<vector<unsigned int> > const &dims) const
{
    return dims[0];
}

}
