#include <config.h>
#include <sarray/SArray.h>
#include "DPar.h"

#include <cmath>
#include <cfloat>

#include <Rmath.h>

using std::vector;

static double ALPHA(vector<SArray const *> const &par)
{
    return *par[0]->value();
}

static double C(vector<SArray const *> const &par)
{
    return *par[1]->value();
}

DPar::DPar()
  : DistReal("dpar", 2, DIST_SPECIAL, true)
{}

DPar::~DPar()
{}

bool DPar::checkParameterValue(vector<SArray const *> const &par) const
{
  return (ALPHA(par) > 0 && C(par) > 0);
}

double 
DPar::d(double x, vector<SArray const *> const &par, bool give_log) const
{
  double alpha = ALPHA(par);
  double c = C(par);

  if (x < c)
    return give_log ? -DBL_MAX : 0;

  if (give_log)
    return log(alpha) + alpha * log(c) - (alpha + 1) * log(x);
  else
    return alpha * exp(alpha * log(c) - (alpha + 1) * log(x));
}

double 
DPar::p(double x, vector<SArray const *> const &par, bool lower, bool give_log)
  const
{
  double alpha = ALPHA(par);
  double c = C(par);
  
  if (x < c)
    return give_log ? -DBL_MAX : 0;
  
  double logq = alpha * log(c/x);
  if (!lower) {
    return give_log ? logq : exp(logq);
  }
  else {
    return give_log ? log(1 - exp(logq)) : 1 - exp(logq);
  }
}

double 
DPar::q(double p, vector<SArray const *> const &par, bool lower, 
	bool log_p) const
{
    if ((log_p  && p > 0) ||                    
	(!log_p && (p < 0 || p > 1)) )          
	return C(par); /* FIXME: should return NaN */
    
    double logp;

    if (lower) {
	if (log_p)
	    logp = p;
	else
	    logp = log(p);
    }
    else {
	if (log_p)
	    logp = log(1 - exp(p)); 
	else
	    logp = log(1 - p);
    }
  
    return exp(log(C(par)) - logp/ALPHA(par));
}

double DPar::r(vector<SArray const *> const &par, RNG *rng) const
{
    return exp(log(C(par)) - rexp(1, rng)/ALPHA(par));
}

double DPar::l(vector<SArray const*> const &par, bool absolute) const
{
  if (absolute) {
    if (par[1]->isFixed()) {
      return C(par);
    }
    else {
      return 0;
    }
  }
  else {
    return C(par);
  }
}

double DPar::u(vector<SArray const*> const &par, bool absolute) const
{
  return DBL_MAX;
}
