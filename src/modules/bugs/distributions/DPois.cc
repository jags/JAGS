#include <config.h>
#include <sarray/SArray.h>
#include "DPois.h"

#include <limits.h>
#include <algorithm>

#include <Rmath.h>

using std::vector;
using std::max;

static inline double LAMBDA(vector<SArray const *> const &par)
{
    return *par[0]->value();
}

DPois::DPois()
  : DistDiscrete("dpois", 1, DIST_POSITIVE, true)
{}

DPois::~DPois()
{}

bool DPois::checkParameterValue(vector<SArray const *> const &par) const
{
    return (LAMBDA(par) > 0);
}

double
DPois::d(double x, vector<SArray const *> const &par, bool give_log) const
{
    return dpois(x, LAMBDA(par), give_log);
}

double
DPois::p(double q, vector<SArray const *> const &par, bool lower, bool give_log)
    const
{
    return ppois(q, LAMBDA(par), lower, give_log);
}

double 
DPois::q(double p, vector<SArray const *> const &par, bool lower, bool log_p)
    const
{
    return qpois(p, LAMBDA(par), lower, log_p);
}

double DPois::r(vector<SArray const *> const &par, RNG *rng) const
{
    return rpois(LAMBDA(par), rng);
}
