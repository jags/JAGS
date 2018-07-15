#include <config.h>
#include "DNegBin.h"
#include <util/nainf.h>

#include <cmath>
#include <vector>
#include <algorithm>

#include <JRmath.h>

using std::vector;
using std::max;
using std::string;

#define PROB(par) (*par[0])
#define SIZE(par) (*par[1])

namespace jags {
namespace bugs {

DNegBin::DNegBin()
    : RScalarDist("dnegbin", 2, DIST_POSITIVE, true)
{}

string DNegBin::alias() const
{
    return "dnbinom";
}

bool 
DNegBin::checkParameterValue (vector<double const *> const &par) const
{
    double p = PROB(par);
    return (SIZE(par) >= 0 && p > 0 && p <= 1);
}

double
DNegBin::d(double x, PDFType,
	   vector<double const *> const &par, bool give_log) 
    const
{
    if (SIZE(par) == 0) {
	if (give_log) {
	    return (x == 0) ? 0 : JAGS_NEGINF;
	}
	else {
	    return (x == 0) ? 1 : 0;
	}
    }
    else {
	return dnbinom(x, SIZE(par), PROB(par), give_log);
    }
}

double
DNegBin::p(double q, vector<double const *> const &par, bool lower, 
	   bool give_log) const
{
    if (SIZE(par) == 0) {
	return give_log ? 0 : 1;
    }
    else {
	return pnbinom(q, SIZE(par), PROB(par), lower, give_log);
    }
}

double 
DNegBin::q(double p, vector<double const *> const &par, bool lower, 
	   bool log_p) const
{
    if (SIZE(par) == 0) {
	return 0;
    }
    else {
	return qnbinom(p, SIZE(par), PROB(par), lower, log_p);
    }
}

double DNegBin::r(vector<double const *> const &par, RNG *rng) const
{
    if (SIZE(par) == 0) {
	return 0;
    }
    else {
	return rnbinom(SIZE(par), PROB(par), rng);
    }
}

double DNegBin::KL(vector<double const *> const &par0,
		   vector<double const *> const &par1) const
{
    double p0 = PROB(par0);
    double r0 = SIZE(par0);
    double p1 = PROB(par1);
    double r1 = SIZE(par1);

    if (fabs(r0 - r1) > 1e-16) {
	//We can't calculat Kullback-Leibler divergence in closed form when
	//r0 and r1 are different
	return JAGS_NA;
    }
    
    return r0 * (log(p0) - log(p1))  + 
	(1 - p0) * r0 * (log(1 - p0) - log(1 - p1)) / p0;
}

}}
