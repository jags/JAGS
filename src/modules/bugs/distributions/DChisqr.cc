#include <config.h>
#include <sarray/SArray.h>
#include "DChisqr.h"

#include <Rmath.h>

using std::vector;

static inline double DF (vector<SArray const *> const &par)
{
    return *par[0]->value();
}

DChisqr::DChisqr()
  : DistReal("dchisqr", 1, DIST_POSITIVE, true)
{}

DChisqr::~DChisqr()
{}


bool 
DChisqr::checkParameterValue(vector<SArray const *> const &par) const
{
  return (DF(par) > 0);
}

double 
DChisqr::d(double x, vector<SArray const *> const &par, bool give_log) const
{
  return dchisq(x, DF(par), give_log);
}

double 
DChisqr::p(double q, vector<SArray const *> const &par, bool lower, bool log_p)
  const
{
  return pchisq(q, DF(par), lower, log_p);
}

double
DChisqr::q(double p, vector<SArray const *> const &par, bool lower, bool log_p)
const
{
  return qchisq(p, DF(par), lower, log_p);
}

double DChisqr::r(vector<SArray const *> const &par, RNG *rng) const
{
    return rchisq(DF(par), rng);
}
