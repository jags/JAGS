#include <config.h>
#include <util/logical.h>
#include <util/nainf.h>
#include <util/dim.h>
#include <module/ModuleError.h>

#include "DSum.h"

#include <cfloat>
#include <cmath>

using std::vector;
using std::fabs;
using std::sqrt;

namespace jags {
namespace bugs {

DSum::DSum()
    : ArrayDist("dsum", 0)
{
}

bool DSum::isDiscreteValued(vector<bool> const &mask) const
{
    return allTrue(mask);
}

double DSum::logDensity(double const *x, PDFType,
			vector<double const *> const &par,
			vector<vector<unsigned long> > const &dims) const
{
    const double tol = sqrt(DBL_EPSILON);
    unsigned long length = product(dims[0]);
    for (unsigned long i = 0; i < length; ++i) {
	double s = x[i];
	for (unsigned long j = 0; j < par.size(); ++j) {
	    s -= par[j][i];
	}
	if (fabs(s) > tol) {
	    return JAGS_NEGINF;
	}
    }
    return 0;
}

void DSum::randomSample(double *x,
			vector<double const *> const &par, 
			vector<vector<unsigned long> > const &dims,
			RNG *) const
{
    unsigned long length = product(dims[0]);
    for (unsigned long i = 0; i < length; ++i) {
	x[i] = 0;
	for (unsigned long j = 0; j < par.size(); ++j) {
	    x[i] += par[j][i];
	}
    }
}

bool DSum::isSupportFixed(vector<bool> const &fixmask) const
{
    return allTrue(fixmask);
}

    bool DSum::fullRank() const
    {
	return false;
    }

bool DSum::checkParameterValue(vector<double const *> const &,
			       vector<vector<unsigned long> > const &) const
{
    return true;
}

bool DSum::checkParameterDim (vector<vector<unsigned long> > const &dims) const
{
    if (dims.empty()) return false;
    if (isFlat(dims[0])) return false;
    for (unsigned long i = 1; i < dims.size(); ++i) {
	if (dims[i] != dims[0])
	    return false;
    }
    return true;
}

bool DSum::checkParameterDiscrete(vector<bool> const &mask) const
{
    for (unsigned long i = 1; i < mask.size(); ++i) {
	if (mask[i] != mask[0])
	    return false;
    }
    return true;
}

void DSum::support(double *lower, double *upper,
		   vector<double const *> const &par,
		   vector<vector<unsigned long> > const &dims) const
{
    unsigned long length = product(dims[0]);
    for (unsigned long i = 0; i < length; ++i) {
	lower[i] = 0;
	for (unsigned long j = 0; j < par.size(); ++j) {
	    lower[i] += par[j][i];
	}
	upper[i] = lower[i];
    }
}

vector<unsigned long> DSum::dim(vector<vector<unsigned long> > const &dims) const
{
    return dims[0];
}

}}
