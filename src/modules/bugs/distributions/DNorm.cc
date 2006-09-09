#include <config.h>
#include <sarray/SArray.h>
#include "DNorm.h"

#include <cmath>

#include <Rmath.h>

using std::vector;

static double MU(vector<SArray const *> const &par)
{
    return *par[0]->value();
}

static double SIGMA(vector<SArray const *> const &par)
{
    return 1/sqrt(*par[1]->value());
}

static double TAU(vector<SArray const *> const &par)
{
    return *par[1]->value();
}


DNorm::DNorm()
  : DistReal("dnorm", 2, DIST_UNBOUNDED, true)
{}

DNorm::~DNorm()
{}

bool DNorm::checkParameterValue(vector<SArray const *> const &par) const
{
    return (TAU(par) > 0);
}

double
DNorm::d(double x, vector<SArray const *> const &par, bool give_log) const
{
    return dnorm(x, MU(par), SIGMA(par), give_log);
}

double
DNorm::p(double q, vector<SArray const *> const &par, bool lower, bool give_log)
  const
{
    return pnorm(q, MU(par), SIGMA(par), lower, give_log);
}

double 
DNorm::q(double p, vector<SArray const *> const &par, bool lower, bool log_p)
  const
{
    return qnorm(p, MU(par), SIGMA(par), lower, log_p);
}

double 
DNorm::r(vector<SArray const *> const &par, RNG *rng) const
{
    return rnorm(MU(par), SIGMA(par), rng);
}

