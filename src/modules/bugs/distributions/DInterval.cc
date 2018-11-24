#include <config.h>
#include "DInterval.h"
#include <util/dim.h>
#include <util/nainf.h>

#include <cfloat>
#include <algorithm>

using std::min;
using std::max;
using std::vector;

#define T(par) (*par[0])
#define CUTPOINTS(par) (par[1])
#define NCUT(lengths) (lengths[1])

static unsigned long value(vector<double const *> const &par, unsigned long ncut)
{
    double t = T(par);
    for (unsigned long i = 0; i < ncut; ++i) {
	if (t <= CUTPOINTS(par)[i])
	    return i;
    }
    return ncut;
}

namespace jags {
namespace bugs {

DInterval::DInterval()
    : VectorDist("dinterval", 2)
{
}

bool DInterval::isDiscreteValued(vector<bool> const &) const
{
    return true;
}

bool DInterval::checkParameterLength(vector<unsigned long> const &lengths) const
{
    return lengths[0] == 1 && lengths[1] >= 1;
}

bool DInterval::checkParameterValue(vector<double const *> const &par,
				    vector<unsigned long> const &lengths) 
    const
{
    for (unsigned long i = 1; i < NCUT(lengths); ++i) {
	if (CUTPOINTS(par)[i] <= CUTPOINTS(par)[i-1])
	    return false;
    }
    return true;
}

double 
DInterval::logDensity(double const *y, PDFType,
		      vector<double const *> const &par,
		      vector<unsigned long> const &lengths) const
{
    if (*y < 0)
	return JAGS_NEGINF;
    
    unsigned long x = static_cast<unsigned long>(*y);
    if (x > NCUT(lengths)) {
	return JAGS_NEGINF;
    }
    else {
	double t = T(par);
	if (x > 0 && t <= CUTPOINTS(par)[x-1])
	    return JAGS_NEGINF;
	else if (x < NCUT(lengths) && t > CUTPOINTS(par)[x])
	    return JAGS_NEGINF;
	else
	    return 0;
    }
}

void DInterval::randomSample(double  *x,
			     vector<double const *> const &par,
			     vector<unsigned long> const &lengths, RNG *) const
{
    /* 
       The random sample from DInterval is not random at all,
       but deterministic.
    */
    *x = static_cast<double>(value(par, NCUT(lengths)));
}

    bool DInterval::fullRank() const
    {
	return false;
    }

void DInterval::support(double *lower, double *upper,
			vector<double const *> const &par,
			vector<unsigned long> const &lengths) const
{
    unsigned long y = value(par, NCUT(lengths));    
    *lower = y;
    *upper = y;
}


bool DInterval::isSupportFixed(vector<bool> const &fixmask) const
{
    return fixmask[0] && fixmask[1];
}

unsigned long DInterval::length(vector<unsigned long> const &) const
{
    return 1;
}

double DInterval::KL(vector<double const *> const &par0,
		     vector<double const *> const &par1,
		     vector<unsigned long> const &lengths) const
{
    unsigned long y0 = value(par0, NCUT(lengths));
    unsigned long y1 = value(par1, NCUT(lengths));
    return (y0 == y1) ? 0 : JAGS_POSINF;
}

}}
