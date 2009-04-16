/* 
 * Adapted from R code published in conjunction with:
 * 
 * Liao, J.G. And Rosen, O. (2001) Fast and Stable Algorithms for
 * Computing and Sampling from the Noncentral Hypergeometric
 * Distribution.  The American Statistician 55, 366-369.
 */

#include <config.h>
#include "DHyper.h"
#include <rng/RNG.h>
#include <util/nainf.h>

#include <algorithm>
#include <stdexcept>

#include <JRmath.h>

using std::vector;
using std::max;
using std::min;
using std::logic_error;

DHyper::DHyper()
    : DistScalarRmath("dhyper", 4, DIST_SPECIAL, false, true)
{}

static void
getParameters(int &n1, int &n2, int &m1, double &psi,
	      vector<double const *> const &parameters)
{
    n1 = static_cast<int>(*parameters[0]);
    n2 = static_cast<int>(*parameters[1]);
    m1 = static_cast<int>(*parameters[2]);
    psi = *parameters[3];
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

bool DHyper::checkParameterValue (vector<double const *> const &parameters,
				  vector<vector<unsigned int> > const &dims) 
  const
{
  int n1,n2,m1;
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

static int modeCompute(int n1, int n2, int m1, double psi)
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
  int mode = static_cast<int>(c/q);
  if (mode >= 0 && mode >= m1 - n2 && mode <= n1 && mode <= m1) {
    return mode;
  }
  else {
    return static_cast<int>(q/a);
  }
}

double rfunction(int n1, int n2, int m1, double psi, int i) {
  return psi * (n1 - i + 1) * (m1 - i + 1)/(i * (n2 - m1 + i));
}  

static void density(double *p, int N,
		    int n1, int n2, int m1, double psi)
{
  int ll = max((int) 0, m1 - n2);
  int uu = min(n1, m1);
  if (N != uu - ll + 1) {
    throw logic_error("Length mismatch calculating hypergeometric density");
  }

  int mode = modeCompute(n1, n2, m1, psi);

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

static int 
sampleLowToHigh(int lower_end, double ran, double const *pi, int N)
{
  //fixme: check lower_end > 0 <= N

  for (int i = lower_end; i < N; ++i) {
    if (ran <= pi[i]) 
      return i;
    ran -= pi[i];
  }
  return N - 1;
}

static int 
sampleHighToLow(int upper_end, double ran, double const *pi, int N)
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

static int singleDraw(int n1, int n2, int m1, double psi, 
		       int mode, double const *pi, int N, double ran) 
{
  if (mode == 0) 
    return sampleLowToHigh(0, ran, pi, N);
  if (mode == N - 1)
    return sampleHighToLow(N - 1, ran, pi, N);

  if (ran < pi[mode]) 
    return mode;

  ran -= pi[mode];
  int lower = mode - 1;
  int upper = mode + 1;
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

double DHyper::d(double z, vector<double const *> const &parameters, 
	 bool give_log) const
{
  int n1,n2,m1;
  double psi;
  getParameters(n1, n2, m1, psi, parameters);

  int x = static_cast<int>(z);
  int ll = max((int) 0, m1 - n2);
  int uu = min(n1, m1);

  double den = 0;
  if (x >= ll && x <= uu) {
    int N = uu - ll + 1;
    double *pi = new double[N];
    density(pi, N, n1, n2, m1, psi);
    den = pi[x - ll];
    delete [] pi;
  }

  if (give_log) {
    return den == 0 ? JAGS_NEGINF : log(den);
  }
  else {
    return den;
  }
}

double DHyper::p(double x, vector<double const *> const &parameters, bool lower,
	         bool give_log) const
{
  int n1,n2,m1;
  double psi;
  getParameters(n1, n2, m1, psi, parameters);

  int ll = max((int) 0, m1 - n2);
  int uu = min(n1, m1);

  double sumpi = 0;
  if (x >= ll) {
    if (x >= uu) {
      sumpi = 1;
    }
    else {
      int N = uu - ll + 1;
      double *pi = new double[N];
      density(pi, N, n1, n2, m1, psi);
      for (int i = ll; i <= x; ++i) {
	sumpi += pi[i - ll];
      }
      delete [] pi;
    }
  }

  if (!lower)
    sumpi = max(1 - sumpi, 0.0);
  
  if (give_log) {
    return sumpi == 0 ? JAGS_NEGINF : log(sumpi);
  }
  else {
    return sumpi;
  }
}

double DHyper::q(double p, vector<double const *> const &parameters, bool lower,
               bool log_p) const
{
  int n1,n2,m1;
  double psi;
  getParameters(n1, n2, m1, psi, parameters);
  
  int ll = max((int) 0, m1 - n2);
  int uu = min(n1, m1);
  
  int N = uu - ll + 1;
  double *pi = new double[N];
  density(pi, N, n1, n2, m1, psi);

  if (log_p)
    p = exp(p);
  if (!lower)
    p = 1 - p;

  double sumpi = 0;
  for (int i = ll; i < uu; ++i) {
    sumpi += pi[i - ll];
    if (sumpi >= p) {
      delete [] pi;
      return i;
    }
  }
  delete [] pi;
  return uu;
}

double DHyper::r(vector<double const *> const &parameters, RNG *rng) const
{
  int n1,n2,m1;
  double psi;
  getParameters(n1, n2, m1, psi, parameters);

  int mode = modeCompute(n1, n2, m1, psi);
  int N = max((int) 0, m1 - n2) - min(n1, m1) + 1;
  double *pi = new double[N];
  int y =  singleDraw(n1, n2, m1, psi, mode, pi, N, rng->uniform());
  delete [] pi;
  return y;
}

double DHyper::l(std::vector<double const *> const &parameters) const
{
    int n1,n2,m1;
    double psi;
    getParameters(n1, n2, m1, psi, parameters);
  
    return max((int) 0, m1 - n2);
}

double DHyper::u(std::vector<double const *> const &parameters) const
{
  int n1,n2,m1;
  double psi;
  getParameters(n1, n2, m1, psi, parameters);
  
  return min(n1, m1);
}

bool DHyper::isSupportFixed(vector<bool> const &fixmask) const
{
  return fixmask[0] && fixmask[1] && fixmask[2]; //Margins fixed
}
