#include <config.h>
#include <distribution/DistDiscrete.h>
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
using std::runtime_error;
using std::string;
using std::logic_error;

DistDiscrete::DistDiscrete(string const &name, unsigned int npar, 
			   Support support, bool canbound)
  : DistScalar(name, npar, support, canbound, true)
{
}

double 
DistDiscrete::logLikelihood(SArray const &x,
			    vector<SArray const *> const &parameters) const
{
  bool lb = lowerBound(parameters);
  bool ub = upperBound(parameters);

  if (!lb && !ub) {
    return d(*x.value(), parameters, true);
  }
  else {
    double y = *x.value();
    double loglik =  d(y, parameters, true);
    if (lb && !ub) {
      double ll = lowerSupport(0, parameters);
      if (y < ll) {
	return -DBL_MAX;
      }
      loglik -= p(ll-1, parameters, false, true);
    }
    else if (!lb && ub) {
      double ul = upperSupport(0, parameters);
      if (y > ul) {
	return -DBL_MAX;
      }
      loglik -= p(ul, parameters, true, true);
    }
    else {
      double ll = lowerSupport(0, parameters);
      double ul = upperSupport(0, parameters);
      if (y < ll || y > ul) {
	return -DBL_MAX;
      }
      if (ll > ul) {
	string msg("Inconsistent bounds for distribution ");
	throw runtime_error(msg + name());
      }
      loglik -= log(p(ul, parameters, true, false) -
		    p(ll - 1, parameters, true, false));
    }
    return loglik;
  }
}

void DistDiscrete::randomSample(SArray &x,
				vector<SArray const *> const &parameters,
				RNG *rng) const
{
  if (x.length() != 1)
    throw length_error("length mismatch in DistDiscrete::randomSample");

  bool bb = lowerBound(parameters) != 0;
  bool ba = upperBound(parameters) != 0;
  
  double y;
  if (ba || bb) {
    double plower = p(l(parameters, false) - 1, parameters, true, false);
    double pupper = p(u(parameters, false), parameters, true, false);
    double u = plower + rng->uniform() * (pupper - plower);
    y = q(u, parameters, true, false);
  }
  else {
    y = r(parameters, rng);
  }
  x.setValue(&y,1);
}
