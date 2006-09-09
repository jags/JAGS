#include <config.h>
#include <rng/RNG.h>
#include <sarray/SArray.h>
#include "DCat.h"

#include <cfloat>
#include <algorithm>
#include <cmath>

#include <Rmath.h>

using std::min;
using std::max;
using std::vector;

DCat::DCat()
  : DistDiscrete("dcat", 1, DIST_SPECIAL, true) 
{}

DCat::~DCat()
{}

bool DCat::checkParameterDim(vector<Index> const &dims) const
{
  return (dims[0].size() == 1) && (dims[0][0] > 1);
}

bool DCat::checkParameterValue(vector<SArray const *> const &parameters) const
{
  double const *prob = parameters[0]->value();
  long length = parameters[0]->length();

  double sump = 0.0;
  for (long i = 0; i < length; i++) {
    if (prob[i] < 0.0 || prob[i] > 1.0) {
      return false; // Invalid probability
    }
    sump += prob[i];
  }
  if (fabs(sump - 1.0) > 16*DBL_EPSILON)
    return false; // Probabilities do not sum to 1 

  return true;
}

double DCat::l(vector<SArray const *> const &parameters, bool absolute) const
{
   return 1;
}

double DCat::u(vector<SArray const *> const &parameters, bool absolute) const
{
   return parameters[0]->length();
}

double DCat::d(double x, std::vector<SArray const *> const &parameters, 
	       bool give_log) const
{
  double const *prob = parameters[0]->value();
  long length = parameters[0]->length();
  double value = 0;

  long y = static_cast<long>(x);
  if (y >= 1 && y <= length) {
    value =  prob[y-1];
  }

  if (give_log) {
    return value == 0 ? -DBL_MAX : log(value);
  }
  else {
    return value;
  }
}

double DCat::p(double x, std::vector<SArray const *> const &parameters,
               bool lower, bool give_log) const
{
  double const *prob = parameters[0]->value();

  double sump = 0;  
  if (x > parameters[0]->length()) {
    sump = 1;
  }
  else {
    for (long i = 1; i <= x; ++i) {
      sump += prob[i-1];
    }
  }

  if (!lower) {
    sump = max(1 - sump, 0.0);
  }

  if (give_log) {
    return sump == 0 ? -DBL_MAX : log(sump);
  }
  else {
    return sump;
  }
}

double DCat::q(double p, std::vector<SArray const *> const &parameters,
	     bool lower, bool log_p) const
{
  if (log_p)
    p = exp(p);
  if (!lower)
    p = 1 - p;

  long length = parameters[0]->length();
  double const *prob = parameters[0]->value();

  double sump = 0;
  for (long i = 1; i < length; ++i) {
    sump += prob[i-1];
    if (sump > p + 16*DBL_EPSILON)
      return i;
  }
  return length;
}

double DCat::r(std::vector<SArray const *> const &parameters, RNG *rng) const
{
  return q(rng->uniform(), parameters, true, false);
}

