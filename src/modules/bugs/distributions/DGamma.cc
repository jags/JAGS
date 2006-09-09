#include <config.h>
#include <sarray/SArray.h>
#include "DGamma.h"

#include <Rmath.h>

using std::vector;

static inline double SHAPE(vector<SArray const *> const &par)
{
    return *par[0]->value();
}

static inline double SCALE(vector<SArray const *> const &par)
{
    return 1/(*par[1]->value());
}

static inline double RATE(vector<SArray const *> const &par)
{
    return *par[1]->value();
}

DGamma::DGamma()
  : DistReal("dgamma", 2, DIST_POSITIVE, true)
{}

DGamma::~DGamma()
{}

bool 
DGamma::checkParameterValue(vector<SArray const *> const &par) const
{
    return (SHAPE(par) > 0 && RATE(par) > 0);
}

double
DGamma::d(double x, vector<SArray const *> const &par, bool give_log) 
    const
{
    return dgamma(x, SHAPE(par), SCALE(par), give_log);
}

double
DGamma::p(double q, vector<SArray const *> const &par, bool lower,
	  bool give_log) const
{
    return pgamma(q, SHAPE(par), SCALE(par), lower, give_log);
}

double 
DGamma::q(double p, vector<SArray const *> const &par, bool lower, 
	  bool log_p) const
{
    return qgamma(p, SHAPE(par), SCALE(par), lower, log_p);
}

double DGamma::r(vector<SArray const *> const &par, RNG *rng) const
{
    return rgamma(SHAPE(par), SCALE(par), rng);
}

void DGamma::typicalValue(SArray &x, vector<SArray const *> const &par) const
{
    double mean = SHAPE(par) * SCALE(par);
    x.setValue(&mean, 1);
}

