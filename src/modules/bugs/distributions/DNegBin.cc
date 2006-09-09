#include <config.h>
#include <sarray/SArray.h>
#include "DNegBin.h"

#include <cmath>
#include <vector>
#include <algorithm>

#include <Rmath.h>

using std::vector;
using std::max;

static inline double PROB(vector<SArray const *> const &par)
{
    return *par[0]->value();
}

static inline double SIZE(vector<SArray const *> const &par)
{
    return *par[1]->value();
}

DNegBin::DNegBin()
    : DistDiscrete("dnegbin", 2, DIST_POSITIVE, true)
{}

DNegBin::~DNegBin()
{}

bool DNegBin::checkParameterDiscrete (vector<bool> const &mask) const
{
    return mask[1]; // r is discrete-valued
}

bool 
DNegBin::checkParameterValue(vector<SArray const *> const &par) const
{
    double p = PROB(par);
    return (SIZE(par) > 0 && p > 0 && p < 1);
}

double
DNegBin::d(double x, vector<SArray const *> const &par, bool give_log) 
    const
{
  return dnbinom(x, SIZE(par), PROB(par), give_log);
}

double
DNegBin::p(double q, vector<SArray const *> const &par, bool lower, 
	   bool give_log) const
{
    return pnbinom(q, SIZE(par), PROB(par), lower, give_log);
}

double 
DNegBin::q(double p, vector<SArray const *> const &par, bool lower, 
	   bool log_p) const
{
    return qnbinom(p, SIZE(par), PROB(par), lower, log_p);
}

double DNegBin::r(vector<SArray const *> const &par, RNG *rng) const
{
    return rnbinom(SIZE(par), PROB(par), rng);
}

