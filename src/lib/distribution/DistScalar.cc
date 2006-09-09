#include <config.h>
#include <distribution/DistScalar.h>
#include <sarray/SArray.h>
#include <sarray/util.h>

#include <stdexcept>
#include <cfloat>
#include <cmath>
#include <algorithm>

using std::max;
using std::min;
using std::string;
using std::vector;
using std::length_error;
using std::logic_error;

DistScalar::DistScalar(string const &name, unsigned int npar, Support support, 
		       bool canbound, bool discrete)
  
  : Distribution(name, npar, canbound, discrete), 
    _support(support)
{
}

double 
DistScalar::lowerSupport(unsigned long i,
			 vector<SArray const *> const &parameters) const
{
  if (i != 0) {
    throw logic_error("Invalid index in Distreal::lowerSupport");
  }

  SArray const *lb = lowerBound(parameters);
  if (lb) {
    return max(l(parameters, false), *lb->value());
  }
  else {
    return l(parameters, false);
  }
}

double
DistScalar::upperSupport(unsigned long i,
			 vector<SArray const *> const &parameters) const
{
  if (i != 0) {
    throw logic_error("Invalid index in Distreal::upperSupport");
  }

  SArray const *ub = upperBound(parameters);
  if (ub) {
    return min(u(parameters, false), *ub->value());
  }
  else {
    return u(parameters, false);
  }
}

Index DistScalar::dim(vector<Index> const &parameters) const
{
    static const Index ScalarIndex(1);
    return ScalarIndex;
}

double DistScalar::l(vector<SArray const *> const &parameters,
		     bool absolute) const
{
  switch(_support) {
  case DIST_UNBOUNDED:
      return -DBL_MAX;
      break;
  case DIST_POSITIVE: case DIST_PROPORTION:
      return 0;
      break;
  case DIST_SPECIAL:
      //You must overload this function 
      throw logic_error("Cannot call DistScalar::l for special distribution");
  }
  return 0; //Wall
      
}

double DistScalar::u(vector<SArray const *> const &parameters,
		     bool absolute) const
{
  switch(_support) {
  case DIST_UNBOUNDED: case DIST_POSITIVE:
    return DBL_MAX;
    break;
  case DIST_PROPORTION:
    return 1;
    break;
  case DIST_SPECIAL:
    //You must overload this function 
    throw logic_error("Cannot call DistScalar::u for special distribution");
  }
  return 0; //Wall
}

void
DistScalar::typicalValue(SArray &x,
                         std::vector<SArray const *> const &parameters) const
{

    SArray const *bb = lowerBound(parameters);
    SArray const *ba = upperBound(parameters);

    double y;
    if (!ba && !bb) {
      y = q(0.5, parameters, true, false);
    }
    else if (bb && ba) {

      //Find lower bound, upper bound and median
      double lower = lowerSupport(0, parameters);
      double plower = p(lower, parameters, false, false);
      double upper = upperSupport(0, parameters);
      double pupper = p(upper, parameters, false, false);
      double px = (plower + pupper)/2;
      double med = q(px, parameters, true, false);
      
      //Calculate the log densities
      double dlower = d(lower, parameters, true);
      double dupper = d(upper, parameters, true);
      double dmed = d(med, parameters, true);
      
      //Pick the point with the highest density
      if (dmed > dlower && dmed > dupper) {
	y = med;
      }
      else if (dupper > dlower) {
	y = upper;
      }
      else {
	y = lower;
      }
    }
    else if (bb) {
      double lower = lowerSupport(0, parameters);
      double plower = p(lower, parameters, false, false);
      double med = q(plower/2, parameters, false, false);
      
      //Calculate the log densities
      double dlower = d(lower, parameters, true);
      double dmed = d(med, parameters, true);
      
      //Pick the point with the highest density
      if (dmed > dlower) {
	y = med;
      }
      else {
	y = lower;
      }
    }
    else if (ba) {
      double upper = upperSupport(0, parameters);
      double pupper = p(upper, parameters, true, false);
      double med = q(pupper/2, parameters, true, false);
      
      //Calculate the densities
      double dupper = d(upper, parameters, true);
      double dmed = d(med, parameters, true);
      
      //Pick the point with the highest density
      if (dmed > dupper) {
	y = med;
      }
      else {
	y = upper;
      }
    }

    x.setValue(&y,1);
    //x.setValue(q(0.5, parameters, true, false),0);
}

bool DistScalar::checkParameterDim (vector<Index> const &dims) const
{
  return count_if(dims.begin(), dims.end(), isScalar) == dims.size();
}

