#include <config.h>
#include <sarray/SArray.h>
#include "DBeta.h"

#include <Rmath.h>

using std::vector;

static inline double A(vector<SArray const *> const &parameters)
{
  return *parameters[0]->value();
}

static inline double B(vector<SArray const *> const &parameters)
{
  return *parameters[1]->value();
}

DBeta::DBeta()
  : DistReal("dbeta", 2, DIST_PROPORTION, true)
{}

DBeta::~DBeta()
{}

bool DBeta::checkParameterValue(vector<SArray const *> const &par) const
{
    return (A(par) > 0.0 && B(par) > 0.0);
}

double 
DBeta::d(double x, vector<SArray const *> const &par, bool give_log) const
{
    return dbeta(x, A(par), B(par), give_log);
}

double 
DBeta::p(double q, vector<SArray const *> const &par, bool lower, bool log_p) 
  const
{
    return pbeta(q, A(par), B(par), lower, log_p);
}

double 
DBeta::q(double p, vector<SArray const *> const &par, bool lower, bool log_p) 
  const
{
    return qbeta(p, A(par), B(par), lower, log_p);
}

double DBeta::r(vector<SArray const *> const &par, RNG *rng) const
{
    return rbeta(A(par), B(par), rng);
}
