#include <config.h>
#include "DUnif.h"

#include <cmath>
#include <rng/RNG.h>
#include <util/nainf.h>

using std::vector;
using std::log;

#define LOWER(par) (*par[0])
#define UPPER(par) (*par[1])

namespace jags {
namespace bugs {

DUnif::DUnif()
    : ScalarDist("dunif", 2, DIST_SPECIAL)
{}


bool  DUnif::checkParameterValue (vector<double const *> const &par) const
{
    return (LOWER(par) < UPPER(par));
}

double DUnif::logDensity(double x, PDFType type,
			 vector<double const *> const &par,
			 double const *, double const *) const
{
    if (x < LOWER(par) || x > UPPER(par))
	return JAGS_NEGINF;
    else if (type == PDF_PRIOR) 
	return 0;
    else
	return - log(UPPER(par) - LOWER(par));
}

double DUnif::randomSample(vector<double const *> const &par, 
			   double const *, double const *,
			   RNG *rng) const
{
    return LOWER(par) + rng->uniform() * (UPPER(par) - LOWER(par));
}

double DUnif::l(vector<double const*> const &par) const
{
    return LOWER(par);
}

double DUnif::u(vector<double const*> const &par) const
{
    return UPPER(par);
}

bool DUnif::isSupportFixed(vector<bool> const &fixmask) const
{
    return fixmask[0] && fixmask[1]; //Lower and upper bounds fixed
}

double DUnif::KL(vector<double const *> const &par0,
		 vector<double const *> const &par1) const
{
    if (LOWER(par0) < LOWER(par1))
	return JAGS_POSINF;
    if (UPPER(par0) > UPPER(par1))
	return JAGS_POSINF;

    return log(UPPER(par1) - LOWER(par1)) - log(UPPER(par0) - LOWER(par0));
}

}}
