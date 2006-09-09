#include <config.h>
#include <sarray/SArray.h>
#include "DInterval.h"

#include <cfloat>
#include <algorithm>
#include <stdexcept>

using std::min;
using std::max;
using std::logic_error;
using std::vector;

static unsigned long value(vector<SArray const *> const &parameters)
{
    double t = *parameters[0]->value();
    double const *cutpoints = parameters[1]->value();
    unsigned long length = parameters[1]->length();

    for (unsigned long i = 0; i < length; ++i) {
	if (t <= cutpoints[i])
	    return i;
    }
    return length;
}

DInterval::DInterval()
  : Distribution("dinterval", 2, false, true)
{
}

DInterval::~DInterval()
{}

bool 
DInterval::checkParameterDim(vector<Index> const &dims) const
{
  return dims[0].size() == 1 && dims[0][0] == 1 && dims[1].size() == 1;
}

Index DInterval::dim(vector<Index> const &dims) const
{
    return Index(1);
}

bool DInterval::checkParameterValue(vector<SArray const *> const &parameters) const
{
  unsigned long length = parameters[1]->length();
  double const *cutpoints = parameters[1]->value();
  for (unsigned long i = 1; i < length; ++i) {
    if (cutpoints[i] < cutpoints[i-1])
      return false;
  }
  return true;
}

double 
DInterval::logLikelihood(SArray const &y, 
			 vector<SArray const *> const &parameters) const
{
  if (y.value() < 0)
    return -DBL_MAX;

  unsigned long x = static_cast<unsigned long>(*y.value());
  if (x > parameters[1]->length()) {
    return -DBL_MAX;
  }
  else {
    double t = *parameters[0]->value();
    double const *cutpoints = parameters[1]->value();
    
    if (x > 0 && t <= cutpoints[x-1])
      return -DBL_MAX;
    else if (x < parameters[1]->length() && t > cutpoints[x])
      return -DBL_MAX;
    else
      return 0;
  }
}

void DInterval::randomSample(SArray  &x,
			     vector<SArray const *> const &parameters,
			     RNG *rng) const
{
    /* 
       The random sample from DInterval is not random at all,
       but deterministic.
    */
    x.setValue(static_cast<double>(value(parameters)), 0);
}

unsigned long 
DInterval::df(std::vector<SArray const *> const &parameters) const
{
    return 0;
}

double 
DInterval::lowerSupport(unsigned long i,
			std::vector<SArray const *> const &parameters) const
{
  if (i != 0)
    throw logic_error("Invalid index in DInterval::lowerSupport");

  return value(parameters);
}

double 
DInterval::upperSupport(unsigned long i,
			std::vector<SArray const *> const &parameters) const
{
  if (i != 0)
    throw logic_error("Invalid index in DInteval::upperSupport");
    
  return value(parameters);
}

void
DInterval::typicalValue(SArray &x, 
			std::vector<SArray const *> const &parameters) const
{
     x.setValue(static_cast<double>(value(parameters)), 0);
}
