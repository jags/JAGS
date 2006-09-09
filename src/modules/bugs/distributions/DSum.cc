#include <config.h>
#include <sarray/util.h>
#include <sarray/SArray.h>
#include "DSum.h"

#include <cfloat>
#include <cmath>
#include <stdexcept>

using std::vector;
using std::fabs;
using std::runtime_error;
using std::logic_error;

static inline double SUM(vector<SArray const *> const &par)
{
    return *par[0]->value() + *par[1]->value();
}

DSum::DSum()
  : Distribution("dsum", 2, false, true)
{
}

DSum::~DSum()
{}

bool DSum::checkParameterValue(vector<SArray const *> const &par) const
{
  return true;
}

bool DSum::checkParameterDim(vector<Index> const &dims) const
{
    return isScalar(dims[0]) && isScalar(dims[1]);
}

bool DSum::checkParameterDiscrete(vector<bool> const &mask) const
{
    return allTrue(mask);
}

Index DSum::dim(vector<Index> const &dims) const
{
    return Index(1);
}

double DSum::logLikelihood(SArray const &x, 
			   vector<SArray const *> const &par) const
{
    if (fabs(*x.value() - SUM(par)) > 16*DBL_EPSILON) {
	// If this happens by accident, you have no chance of getting it right
	throw runtime_error("Inconsistent arguments for dsum");
    }
    
    return 0;
}

void
DSum::randomSample(SArray &x, vector<SArray const *> const &par,
		   RNG *rng) const
{
    /* The random sample from DSum is not random at all, but
       deterministic. */
  
    x.setValue(SUM(par), 0);
}

unsigned long DSum::df(std::vector<SArray const *> const &par) const
{
    return 0;
}

double 
DSum::lowerSupport(unsigned long i,
		   std::vector<SArray const *> const &par) const
{
    if (i != 0)
	throw logic_error("Invalid index in DSum::lowerSupport");
    
    return SUM(par);
}


double 
DSum::upperSupport(unsigned long i,
		   std::vector<SArray const *> const &par) const
{
    if (i != 0)
	throw logic_error("Invalid index in DSum::upperSupport");
    
    return SUM(par);
}

void
DSum::typicalValue(SArray &x, 
		   std::vector<SArray const *> const &par) const
{
  x.setValue(SUM(par), 0);
}
