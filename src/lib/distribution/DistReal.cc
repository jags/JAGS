#include <config.h>
#include <distribution/DistReal.h>
#include <rng/RNG.h>
#include <sarray/SArray.h>

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


DistReal::DistReal(string const &name, unsigned int npar, Support support, 
		   bool canbound)
  : DistScalar(name, npar, support, canbound, false)
{
}

double DistReal::logLikelihood(SArray const &x,
			       vector<SArray const *> const &parameters) const
{
  SArray const *lb = lowerBound(parameters);
  SArray const *ub = upperBound(parameters);

  double y = x.value()[0];
  double loglik = d(y, parameters, true);
  if (lb && ub) {
    double ll = lowerSupport(0, parameters);
    if (y < ll)
      return -DBL_MAX;
    double uu = upperSupport(0, parameters); 
    if (y > uu) {
      return -DBL_MAX;
    }
    loglik -= log(p(uu, parameters, true, false) - 
		  p(ll, parameters, true, false));
  }
  else if (lb) {
    double ll = lowerSupport(0, parameters);
    if (y < ll)
      return -DBL_MAX;
    loglik -= p(ll, parameters, false, true);
  }
  else if (ub) {
    double uu = upperSupport(0, parameters); 
    if (y > uu) {
      return -DBL_MAX;
    }
    loglik -= p(uu, parameters, true, true);
  }
  return loglik;
}

void 
DistReal::randomSample(SArray &x, vector<SArray const *> const &parameters,
		       RNG *rng) const
{

    SArray const *bb = lowerBound(parameters);
    SArray const *ba = upperBound(parameters);

    double y;
    if (!ba && !bb) {
      y = r(parameters, rng);
    }
    else if (bb && ba) {
	double lower = lowerSupport(0, parameters);
	double plower = p(lower, parameters, false, false);
	double upper = upperSupport(0, parameters);
	double pupper = p(upper, parameters, false, false);
	double px = plower + rng->uniform() * (pupper - plower);
	y = q(px, parameters, true, false);
    }
    else if (bb) {
	double lower = lowerSupport(0, parameters);
	double plower = p(lower, parameters, false, false);
	double px = plower * rng->uniform();
//	y = q(px, parameters, false, true);
  	y = q(px, parameters, false, false);
    }
    else if (ba) {
	double upper = upperSupport(0, parameters);
	double pupper = p(upper, parameters, true, false);
	double px = pupper * rng->uniform();
	y = q(px, parameters, true, false);
    }
    x.setValue(&y,1);
}
