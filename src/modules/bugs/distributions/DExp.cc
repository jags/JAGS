#include <config.h>
#include "DExp.h"

#include <cmath>
#include <algorithm>

#include <JRmath.h>

using std::max;
using std::vector;

namespace jags {
namespace bugs {

DExp::DExp()
    : RScalarDist("dexp", 1, DIST_POSITIVE)
{}

/* R functions are parameterized in terms of scale = 1/rate */
static inline double SCALE(vector<double const *> const &par)
{
    return 1 /(*par[0]);
}

bool DExp::checkParameterValue (vector<double const *> const &par) const
{
    return (*par[0] > 0);
}

double
DExp::d(double x, PDFType,
	vector<double const *> const &par, bool give_log)  const
{
    return dexp(x, SCALE(par), give_log);
}

double 
DExp::p(double q, vector<double const *> const &par, bool lower, 
	bool log_p) const
{
    return pexp(q, SCALE(par), lower, log_p);
}

double 
DExp::q(double p, vector<double const *> const &par, bool lower,  
	bool log_p) const
{
    return qexp(p, SCALE(par), lower, log_p);
}

double DExp::r(vector<double const *> const &par, RNG *rng) const
{
    return rexp(SCALE(par), rng);
}


double DExp::KL(vector<double const *> const &par0,
		vector<double const *> const &par1) const
{
    double r = (*par1[0]) / (*par0[0]);
    return r - 1 - log(r);
}

}}
