#include <config.h>
#include <sarray/SArray.h>
#include "DBin.h"

#include <algorithm>

#include <Rmath.h>

using std::vector;
using std::max;
using std::min;

/* BUGS parameterization is in opposite order to R parameterization */

static inline double SIZE(vector<SArray const *> const &par)
{
    return *par[1]->value();
}

static inline double PROB(vector<SArray const *> const &par)
{
    return *par[0]->value();
}

DBin::DBin()
  : DistDiscrete("dbin", 2, DIST_SPECIAL, false)
{}

DBin::~DBin()
{}

bool 
DBin::checkParameterDiscrete (vector<bool> const &mask) const
{
    return mask[1] == true;
}

bool DBin::checkParameterValue(vector<SArray const *> const &par) const
{
    return (SIZE(par) >= 1 && PROB(par) >= 0.0 && PROB(par) <= 1.0);
}

double DBin::d(double x, std::vector<SArray const *> const &par, 
	       bool give_log) const
{
    return dbinom(x, SIZE(par), PROB(par), give_log);
}

double DBin::p(double x, std::vector<SArray const *> const &par, 
	       bool lower, bool give_log) const
{
    return pbinom(x, SIZE(par), PROB(par), lower, give_log);
}

double DBin::q(double p, std::vector<SArray const *> const &par, 
	       bool lower, bool log_p) const
{
    return pbinom(p, SIZE(par), PROB(par), lower, log_p);
}

double DBin::r(std::vector<SArray const *> const &par, RNG *rng) const
{
    return rbinom(SIZE(par), PROB(par), rng);
}

double 
DBin::l(vector<SArray const *> const &par, bool absolute) const
{
    return 0;
}
	

double 
DBin::u(vector<SArray const *> const &par, bool absolute) const
{
  if (absolute) {
    /* Largest possible value that u can take for any parameters */
    if (par[1]->isFixed()) {
      return SIZE(par);
    }
    else {
      return DBL_MAX;
    }
  }
  else {
    /* Largest possible value that u can take for any parameters */
    return SIZE(par);
  }
}
