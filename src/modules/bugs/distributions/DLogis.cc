#include <config.h>
#include <sarray/SArray.h>
#include "DLogis.h"

#include <Rmath.h>

using std::vector;

static inline double MU(vector<SArray const *> const &parameters)
{
    return *parameters[0]->value();
}

static inline double SCALE(vector<SArray const *> const &parameters)
{
    return 1/(*parameters[1]->value());
}

static inline double TAU(vector<SArray const *> const &parameters)
{
    return *parameters[1]->value();
}

DLogis::DLogis()
    : DistReal("dlogis", 2, DIST_UNBOUNDED, true)
{}

DLogis::~DLogis()
{}

bool DLogis::checkParameterValue (vector<SArray const *> const &par) const
{
    return (TAU(par) > 0);
}

double 
DLogis::d(double x, vector<SArray const *> const &par, bool give_log) const
{
    return dlogis(x, MU(par), SCALE(par), give_log);
}

double 
DLogis::p(double q, vector<SArray const *> const &par, bool lower, 
	  bool give_log) const
{
    return plogis(q, MU(par), SCALE(par), lower, give_log);
}

double 
DLogis::q(double p, vector<SArray const *> const &par, bool lower, bool log_p)
    const
{
    return qlogis(p, MU(par), SCALE(par), lower, log_p);
}

double 
DLogis::r(vector<SArray const *> const &par, RNG *rng) const
{
    double ans = rlogis(MU(par), SCALE(par), rng);
    return ans;
}
