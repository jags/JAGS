#include <config.h>
#include <sarray/SArray.h>
#include "DExp.h"

#include <cmath>
#include <algorithm>

#include <Rmath.h>

using std::max;
using std::vector;

DExp::DExp()
  : DistReal("dexp", 1, DIST_POSITIVE, true)
{}

DExp::~DExp()
{}

/* R functions are parameterized in terms of scale = 1/rate */
static inline double SCALE(vector<SArray const *> const &par)
{
    return 1 /(*par[0]->value());
}

bool DExp::checkParameterValue(vector<SArray const *> const &par) const
{
  return (*par[0]->value() > 0);
}

double
DExp::d(double x, vector<SArray const *> const &par, bool give_log) 
    const
{
    return dexp(x, SCALE(par), give_log);
}

double 
DExp::p(double q, vector<SArray const *> const &par, bool lower, 
	bool log_p) const
{
    return pexp(q, SCALE(par), lower, log_p);
}

double 
DExp::q(double p, vector<SArray const *> const &par, bool lower,  
	bool log_p) const
{
    return qexp(p, SCALE(par), lower, log_p);
}

double DExp::r(vector<SArray const *> const &par, RNG *rng) const
{
    return rexp(SCALE(par), rng);
}

