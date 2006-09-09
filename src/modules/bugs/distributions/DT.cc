#include <config.h>
#include <sarray/SArray.h>
#include "DT.h"

#include <cmath>

#include <Rmath.h>

using std::vector;

// Mean
static inline double MU (vector<SArray const *> const &par)
{
    return *par[0]->value();
}
// Precision
static inline double TAU (vector<SArray const *> const &par)
{
    return *par[1]->value();
}
// Degrees of freedom 
static inline double DF (vector<SArray const *> const &par)
{
    return *par[2]->value();
}

DT::DT()
  : DistReal("dt", 3, DIST_UNBOUNDED, true)
{}

DT::~DT()
{}

bool DT::checkParameterValue (vector<SArray const *> const &par) const
{
    return (TAU(par) > 0 && DF(par) > 0);
}

double DT::d(double x, vector<SArray const *> const &par, bool give_log) const
{
    x = (x - MU(par)) * sqrt(TAU(par));
    if (give_log) {
	return dt(x, DF(par), 1) + log(TAU(par))/2;
    }
    else {
	return dt(x, DF(par), 0) * sqrt(TAU(par));
    }
}

double DT::p(double x, vector<SArray const *> const &par, bool lower, 
	     bool use_log) const
{
    return pt((x - MU(par)) * sqrt(TAU(par)), DF(par), lower, use_log);
}

double DT::q(double p, vector<SArray const *> const &par, bool lower, 
	     bool log_p) const
{
    return MU(par) + qt(p, DF(par), lower, log_p) / sqrt(TAU(par));
}

double DT::r(vector<SArray const *> const &par, RNG *rng) const
{
    return rt(DF(par), rng) / sqrt(TAU(par)) + MU(par);
}
