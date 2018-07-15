#include <config.h>
#include "DPar.h"
#include <util/nainf.h>
#include <rng/RNG.h>

#include <cmath>
#include <cfloat>

using std::vector;
using std::exp;
using std::log;

#define ALPHA(par) (*par[0])
#define C(par) (*par[1])

namespace jags {
namespace bugs {

DPar::DPar()
    : RScalarDist("dpar", 2, DIST_SPECIAL)
{}

bool DPar::checkParameterValue (vector<double const *> const &par) const
{
  return (ALPHA(par) > 0 && C(par) > 0);
}

double 
DPar::d(double x, PDFType,
vector<double const *> const &par, bool give_log) const
{
  double alpha = ALPHA(par);
  double c = C(par);

  if (x < c)
    return give_log ? JAGS_NEGINF : 0;

  if (give_log)
    return log(alpha) + alpha * log(c) - (alpha + 1) * log(x);
  else
    return alpha * exp(alpha * log(c) - (alpha + 1) * log(x));
}

double 
DPar::p(double x, vector<double const *> const &par, bool lower, bool give_log)
  const
{
  double alpha = ALPHA(par);
  double c = C(par);
  
  if (x < c)
    return give_log ? JAGS_NEGINF : 0;
  
  double logq = alpha * log(c/x);
  if (!lower) {
    return give_log ? logq : exp(logq);
  }
  else {
    return give_log ? log(1 - exp(logq)) : 1 - exp(logq);
  }
}

double 
DPar::q(double p, vector<double const *> const &par, bool lower, 
	bool log_p) const
{
    if ( (log_p  && p > 0) || (!log_p && (p < 0 || p > 1)) )          
	return JAGS_NAN;
    
    double logp;

    if (!lower) {
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

double DPar::r(vector<double const *> const &par, RNG *rng) const
{
    return C(par) * exp(rng->exponential()/ALPHA(par));
}

double DPar::l(vector<double const*> const &par) const
{
    return C(par);
}

double DPar::u(vector<double const*> const &) const
{
  return JAGS_POSINF;
}

bool DPar::isSupportFixed(vector<bool> const &fixmask) const
{
    return fixmask[1]; //Fixed if C is fixed
}

}}
