#include <config.h>
#include <rng/RNG.h>
#include <sarray/SArray.h>
#include "DBern.h"

#include <cmath>
#include <algorithm>

#include <Rmath.h>

using std::vector;
using std::min;
using std::max;

static inline double PROB(vector<SArray const *> const &parameters)
{
  return *parameters[0]->value();
}

DBern::DBern()
  : DistDiscrete("dbern", 1, DIST_PROPORTION, false)
{}

DBern::~DBern()
{}

bool 
DBern::checkParameterValue(vector<SArray const *> const &parameters) const
{
    return  (PROB(parameters) >= 0.0 && PROB(parameters) <= 1.0);
}

double DBern::d(double x, vector<SArray const *> const &parameters,
		bool give_log) const
{
    double d = 0;
    if (x == 1)
	d = PROB(parameters);
    else if (x == 0)
	d = 1 - PROB(parameters);
  
    if (give_log) {
	return d == 0 ? -DBL_MAX : log(d);
    }
    else {
	return d;
    }
}

double DBern::p(double x, vector<SArray const *> const &parameters, bool lower,
		bool give_log) const
{
    return pbinom(x, 1, PROB(parameters), lower, give_log);
}

double DBern::q(double p, std::vector<SArray const *> const &parameters, 
		bool lower, bool log_p) const
{
    return qbinom(p, 1, PROB(parameters), lower, log_p);
}

double DBern::r(std::vector<SArray const *> const &parameters,
		RNG *rng) const
{
  return rng->uniform() < PROB(parameters) ? 1 : 0;
}

