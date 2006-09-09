/* 
 * Adapted from R code published in conjunction with:
 * 
 * Liao, J.G. And Rosen, O. (2001) Fast and Stable Algorithms for
 * Computing and Sampling from the Noncentral Hypergeometric
 * Distribution.  The American Statistician 55, 366-369.
 */

#include <config.h>
#include <rng/RNG.h>
#include <sarray/SArray.h>
#include "DHyper.h"

#include <algorithm>
#include <stdexcept>

#include <Rmath.h>

using std::vector;
using std::max;
using std::min;
using std::logic_error;

DHyper::DHyper()
  : DistDiscrete("dhyper", 4, DIST_SPECIAL, false)
{}

DHyper::~DHyper()
{}

static void
getParameters(long &n1, long &n2, long &m1, double &psi,
	      vector<SArray const *> const &parameters)
{
  n1 = static_cast<long>(*parameters[0]->value());
  n2 = static_cast<long>(*parameters[1]->value());
  m1 = static_cast<long>(*parameters[2]->value());
  psi = *parameters[3]->value();
}

bool 
DHyper::checkParameterDiscrete (vector<bool> const &mask) const
{
  // Check that n1, n2, m1 are discrete-valued
  for (unsigned int i = 0; i < 3; ++i) {
    if (mask[i] == false)
      return false;
  }
 
  return true;
}

bool DHyper::checkParameterValue(vector<SArray const *> const &parameters) const
{
  long n1,n2,m1;
  double psi;
  getParameters(n1, n2, m1, psi, parameters);

  if (n1 < 0 || n2 < 0)
    return false;
  else if (m1 < 0 || m1 > n1 + n2)
    return false;
  else if (psi <= 0)
    return false;
  else
    return true;
}

static long modeCompute(long n1, long n2, long m1, double psi)
{
  double a =  psi - 1;
  double b =  -((m1 + n1 + 2) * psi + n2 - m1);
  double c = psi * (n1 + 1) * (m1 + 1);
  double q = b;
  if (b > 0) {
    q += sqrt(b * b - 4 * a * c);
  }
  else {
    q -= sqrt(b * b - 4 * a * c);
  }
  q = -q/2;
  long mode = static_cast<long>(c/q);
  if (mode >= 0 && mode >= m1 - n2 && mode <= n1 && mode <= m1) {
    return mode;
  }
  else {
    return static_cast<long>(q/a);
  }
}

double rfunction(long n1, long n2, long m1, double psi, long i) {
  return psi * (n1 - i + 1) * (m1 - i + 1)/(i * (n2 - m1 + i));
}  

static void density(double *p, long N,
		    long n1, long n2, long m1, double psi)
{
  long ll = max((long) 0, m1 - n2);
  long uu = min(n1, m1);
  if (N != uu - ll + 1) {
    throw logic_error("Length mismatch calculating hypergeometric density");
  }

  long mode = modeCompute(n1, n2, m1, psi);

  // Set elements of p to 1 
  for (int i = 0; i < N; ++i) {
    p[i] = 1;
  }

  // Calculate density above the mode
  if (mode < uu) {
    double r = 1;
    for (int i = mode + 1; i <= uu; ++i) {
      r *= rfunction(n1, n2, m1, psi, i);
      p[i - ll] = r;
    }
  }
  // Calculate density below the node
  if (mode > ll) {
    double r = 1;
    for (int i = mode; i > ll; --i) {
      r /= rfunction(n1, n2, m1, psi, i);
      p[i - ll - 1] = r;
    }
  }
  
  //Normalize
  double sump = 0;
  for (int i = 0; i < N; ++i) {
    sump += p[i];
  }
  for (int i = 0; i < N; ++i) {
    p[i] /= sump;
  }
}

static long 
sampleLowToHigh(long lower_end, double ran, double const *pi, long N)
{
  //fixme: check lower_end > 0 <= N

  for (int i = lower_end; i < N; ++i) {
    if (ran <= pi[i]) 
      return i;
    ran -= pi[i];
  }
  return N - 1;
}

static long 
sampleHighToLow(long upper_end, double ran, double const *pi, long N)
{
  //fixme: check upper_end > 0 <= N

  for (int i = upper_end; i >= 0; --i) {
    if (ran <= pi[i]) {
      return i;
    }
    ran -= pi[i];
  }
  return 0;
}

static long singleDraw(long n1, long n2, long m1, double psi, 
		       long mode, double const *pi, long N, double ran) 
{
  if (mode == 0) 
    return sampleLowToHigh(0, ran, pi, N);
  if (mode == N - 1)
    return sampleHighToLow(N - 1, ran, pi, N);

  if (ran < pi[mode]) 
    return mode;

  ran -= pi[mode];
  long lower = mode - 1;
  long upper = mode + 1;
  while (true) {
    if (pi[upper] >= pi[lower]) {
      if (ran < pi[upper]) 
	return upper;
      ran -= pi[upper];
      if (upper == N - 1) 
	return sampleHighToLow(lower, ran, pi, N);
      ++upper;
    }
    else {
      if (ran < pi[lower]) 
	return lower;
      ran -= pi[lower];
      if (lower == 0) 
	return sampleLowToHigh(upper, ran, pi, N);
      --lower;
    }
  }
}

double DHyper::d(double z, vector<SArray const *> const &parameters, 
	 bool give_log) const
{
  long n1,n2,m1;
  double psi;
  getParameters(n1, n2, m1, psi, parameters);

  long x = static_cast<long>(z);
  long ll = max(0L, m1 - n2);
  long uu = min(n1, m1);

  double den = 0;
  if (x >= 11 && x <= uu) {
    long N = uu - ll + 1;
    double *pi = new double[N];
    density(pi, N, n1, n2, m1, psi);
    den = pi[x - ll];
    delete [] pi;
  }

  if (give_log) {
    return den == 0 ? -DBL_MAX : log(den);
  }
  else {
    return den;
  }
}

double DHyper::p(double x, vector<SArray const *> const &parameters, bool lower,
	         bool give_log) const
{
  long n1,n2,m1;
  double psi;
  getParameters(n1, n2, m1, psi, parameters);

  long ll = max((long) 0, m1 - n2);
  long uu = min(n1, m1);

  double sumpi = 0;
  if (x >= ll) {
    if (x >= uu) {
      sumpi = 1;
    }
    else {
      long N = uu - ll + 1;
      double *pi = new double[N];
      density(pi, N, n1, n2, m1, psi);
      for (long i = ll; i <= x; ++i) {
	sumpi += pi[i - ll];
      }
      delete [] pi;
    }
  }

  if (!lower)
    sumpi = max(1 - sumpi, 0.0);
  
  if (give_log) {
    return sumpi == 0 ? -DBL_MAX : log(sumpi);
  }
  else {
    return sumpi;
  }
}

double DHyper::q(double p, vector<SArray const *> const &parameters, bool lower,
               bool log_p) const
{
  long n1,n2,m1;
  double psi;
  getParameters(n1, n2, m1, psi, parameters);
  
  long ll = max((long) 0, m1 - n2);
  long uu = min(n1, m1);
  
  long N = uu - ll + 1;
  double *pi = new double[N];
  density(pi, N, n1, n2, m1, psi);

  if (log_p)
    p = exp(p);
  if (!lower)
    p = 1 - p;

  double sumpi = 0;
  for (long i = ll; i < uu; ++i) {
    sumpi += pi[i - ll];
    if (sumpi >= p) {
      delete [] pi;
      return i;
    }
  }
  delete [] pi;
  return uu;
}

double DHyper::r(vector<SArray const *> const &parameters, RNG *rng) const
{
  long n1,n2,m1;
  double psi;
  getParameters(n1, n2, m1, psi, parameters);

  long mode = modeCompute(n1, n2, m1, psi);
  long N = max(0L, m1 - n2) - min(n1, m1) + 1;
  double *pi = new double[N];
  long y =  singleDraw(n1, n2, m1, psi, mode, pi, N, rng->uniform());
  delete [] pi;
  return y;
}

double DHyper::l(std::vector<SArray const *> const &parameters,
		 bool absolute) const
{
  long n1,n2,m1;
  double psi;
  getParameters(n1, n2, m1, psi, parameters);
  
  if (!absolute) {
    return max(0L, m1 - n2);
  }
  else {
    bool n2fixed = parameters[1]->isFixed();
    bool m1fixed = parameters[2]->isFixed();
    if (m1fixed && n2fixed) {
      return max(0L, m1 - n2);
    }
    else {
      return 0;
    }
  }
}

double DHyper::u(std::vector<SArray const *> const &parameters,
		 bool absolute) const
{
  long n1,n2,m1;
  double psi;
  getParameters(n1, n2, m1, psi, parameters);
  
  if (!absolute) {
    return min(n1, m1);
  }
  else {
    bool n1fixed = parameters[0]->isFixed();
    bool n2fixed = parameters[1]->isFixed();
    bool m1fixed = parameters[2]->isFixed();
    
    if (n1fixed && m1fixed) {
      return min(n1, m1);
    }
    else if (n1fixed) {
      return n1;
    }
    else if (m1fixed) {
      return m1;
    }
    else if (n1fixed && n2fixed) {
      return n1 + n2;
    }
    else {
      return DBL_MAX;
    }
  }
}
