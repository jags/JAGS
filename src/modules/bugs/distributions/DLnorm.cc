#include <config.h>
#include <sarray/SArray.h>
#include "DLnorm.h"

#include <cmath>

#include <Rmath.h>

using std::vector;

static inline double MU(vector<SArray const *> const &par)
{
    return *par[0]->value();
}

static inline double TAU(vector<SArray const *> const &par)
{
    return *par[1]->value();
}

static inline double SDLOG(vector<SArray const *> const &par)
{
    return 1/sqrt(*par[1]->value());
}


DLnorm::DLnorm()
  : DistReal("dlnorm", 2, DIST_POSITIVE, true)
{}

DLnorm::~DLnorm()
{}

bool DLnorm::checkParameterValue(vector<SArray const *> const &par) const
{
    return TAU(par) > 0;
}

double 
DLnorm::d(double x, vector<SArray const *> const &par, bool give_log) const
{
    return dlnorm(x, MU(par), SDLOG(par), give_log);
}

double 
DLnorm::p(double q, vector<SArray const *> const &par, bool lower, 
	  bool give_log) const
{
    return plnorm(q, MU(par), SDLOG(par), lower, give_log);
}

double 
DLnorm::q(double p, vector<SArray const *> const &par, bool lower, bool log_p)
    const
{
    return qlnorm(p, MU(par), SDLOG(par), lower, log_p);
}

double 
DLnorm::r(vector<SArray const *> const &par, RNG *rng) const

{
    return rlnorm(MU(par), SDLOG(par), rng);
}
