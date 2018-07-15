#include <config.h>
#include "DBin.h"
#include <util/nainf.h>

#include <algorithm>

#include <JRmath.h>

#include <util/nainf.h>

using std::vector;
using std::max;
using std::min;
using std::string;

/* BUGS parameterization is in opposite order to R parameterization */
#define SIZE(par) (*par[1])
#define PROB(par) (*par[0])

namespace jags {
namespace bugs {

DBin::DBin()
    : RScalarDist("dbin", 2, DIST_SPECIAL, true)
{}

string DBin::alias() const
{
    return "dbinom";
}

bool 
DBin::checkParameterDiscrete (vector<bool> const &mask) const
{
    return mask[1] == true;
}

bool DBin::checkParameterValue (vector<double const *> const &par) const
{
    return (SIZE(par) >= 0 && PROB(par) >= 0.0 && PROB(par) <= 1.0);
}

double DBin::d(double x, PDFType, vector<double const *> const &par, 
	       bool give_log) const
{
    return dbinom(x, SIZE(par), PROB(par), give_log);
}

double DBin::p(double x, vector<double const *> const &par, 
	       bool lower, bool give_log) const
{
    return pbinom(x, SIZE(par), PROB(par), lower, give_log);
}

double DBin::q(double p, vector<double const *> const &par, 
	       bool lower, bool log_p) const
{
    return qbinom(p, SIZE(par), PROB(par), lower, log_p);
}

double DBin::r(vector<double const *> const &par, RNG *rng) const
{
    return rbinom(SIZE(par), PROB(par), rng);
}

double 
DBin::l(vector<double const *> const &) const
{
    return 0;
}
	

double 
DBin::u(vector<double const *> const &par) const
{
    return SIZE(par);
}

bool DBin::isSupportFixed(vector<bool> const &fixmask) const
{
   return fixmask[1]; //SIZE is fixed;
} 

double DBin::KL(vector<double const *> const &par0,
		vector<double const *> const &par1) const
{
    double N0 = SIZE(par0);
    double N1 = SIZE(par1);
    double p0 = PROB(par0);
    double p1 = PROB(par1);

    //FIXME: If N0 < N1 then we can still calculate finite KL
    if (N0 != N1) {
	return JAGS_POSINF;
    }
    else if (p0 == 0) {
	return - N0 * log(1 - p1);
    }
    else if (p0 == 1) {
	return - N0 * log(p1);
    }
    else {
	return (N0 * p0 * (log(p0) - log(p1)) +
		N0 * (1 - p0) * (log(1 - p0) - log(1 - p1)));
    }
}

}}
