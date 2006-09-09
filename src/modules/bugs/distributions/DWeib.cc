#include <config.h>
#include <sarray/SArray.h>
#include "DWeib.h"

#include <vector>
#include <cmath>

#include <Rmath.h>

using std::vector;

static inline double SHAPE (vector<SArray const *> const &par) 
{
    return *par[0]->value();
}

static inline double SCALE (vector<SArray const *> const &par) 
{
    return pow(*par[1]->value(), -1/(*par[0]->value()));
}

DWeib::DWeib()
    : DistReal("dweib", 2, DIST_POSITIVE, true)
{}

DWeib::~DWeib()
{}

bool 
DWeib::checkParameterValue (vector<SArray const *> const &par) const
{
    return (*par[0]->value() > 0 && *par[1]->value() > 0);
}

double 
DWeib::d(double x, vector<SArray const *> const &par, bool give_log)
    const
{
    return dweibull(x, SHAPE(par), SCALE(par), give_log);
}

double
DWeib::p(double q, vector<SArray const *> const &par, bool lower, bool give_log)
  const
{
    return pweibull(q, SHAPE(par), SCALE(par), lower, give_log);
}

double 
DWeib::q(double p, vector<SArray const *> const &par, bool lower, bool log_p)
    const
{
    return qweibull(p, SHAPE(par), SCALE(par), lower, log_p);
}

double DWeib::r(vector<SArray const *> const &par, RNG *rng) const
{
    return rweibull(SHAPE(par), SCALE(par), rng);
}
