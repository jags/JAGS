#include <config.h>
#include "DBern.h"
#include <rng/RNG.h>
#include <util/nainf.h>

#include <cmath>
#include <algorithm>

#include <JRmath.h>

using std::vector;
using std::min;
using std::max;

#define PROB(par) (*par[0])

namespace jags {
namespace bugs {

DBern::DBern()
    : ScalarDist("dbern", 1, DIST_PROPORTION)
{}

bool 
DBern::checkParameterValue (vector<double const *> const &parameters) const
{
    return  (PROB(parameters) >= 0.0 && PROB(parameters) <= 1.0);
}

double DBern::logDensity(double x, PDFType ,
			 vector<double const *> const &parameters,
			 double const *, double const *) const 
{
    double d = 0;
    if (x == 1)
	d = PROB(parameters);
    else if (x == 0)
	d = 1 - PROB(parameters);
    
    return d == 0 ? JAGS_NEGINF : log(d);
}

double DBern::randomSample(vector<double const *> const &parameters, 
			   double const *, double const *,
			   RNG *rng) const
{
    return rng->uniform() < PROB(parameters) ? 1 : 0;
}

bool DBern::canBound() const
{
    return false;
}

bool DBern::isDiscreteValued(vector<bool> const &) const
{
    return true;
}

    double DBern::KL(vector<double const *> const &par0,
		     vector<double const *> const &par1) const
    {
	double p0 = PROB(par0);
	double p1 = PROB(par1);
	
	if (p0 == 0) {
	    return - log(1 - p1);
	}
	else if (p0 == 1) {
	    return - log(p1);
	}
	else {
	    return (p0 * (log(p0) - log(p1)) +
		    (1 - p0) * (log(1 - p0) - log(1 - p1)));
	}
    }

}}
