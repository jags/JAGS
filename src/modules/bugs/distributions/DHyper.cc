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
#include <vector>
#include <stdexcept>
#include <cfloat>

#include <JRmath.h>

using std::vector;
using std::max;
using std::min;
using std::vector;
using std::logic_error;

namespace jags {
namespace bugs {

DHyper::DHyper()
    : RScalarDist("dhyper", 4, DIST_SPECIAL, true)
{}

bool DHyper::canBound() const
{
    return false;
}

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

bool DHyper::checkParameterValue (vector<double const *> const &params) const
{
    int n1,n2,m1;
    double psi;
    getParameters(n1, n2, m1, psi, params);

    if (n1 < 0 || n2 < 0)
	return false;
    else if (m1 < 0 || m1 > n1 + n2)
	return false;
    else if (psi <= 0)
	return false;
    else
	return true;
}

/* 
 * Calculates the mode of the hypergeometric distribution
 *
 * Extending the unnormalized density function to continuous x, We
 * solve the equation p(x) = p(x-1). Using the recurrence relation
 * given in rfunction (see below) this reduces to a quadratic
 * equation. The maximum value for integer values must lie between x-1
 * and x, so to find the mode we take the floor of x.
 */
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

/* 
 * Recurrence relation that allows us to calculate hypergeometric
 * probabilities without running into combinatoric problems:
 *  
 * p(x) = p(x - 1) * rfunction(n1, n2, m1, psi, x)
 */
static inline double rfunction(int n1, int n2, int m1, double psi, int i) 
{
    return psi * (n1 - i + 1) * (m1 - i + 1)/(i * (n2 - m1 + i));
}  

/**
 * Calculates the unnormalized density of x. The mode has reference
 * density value 1. 
 */
static double density_unnormalized(int x, int n1, int n2, int m1, double psi)
{
    int ll = max(0, m1 - n2);
    int uu = min(n1, m1);

    if (x < ll || x > uu) return 0;

    // Density at mode has reference value 1
    int mode = modeCompute(n1, n2, m1, psi);
    double r = 1;

    if (x > mode) {
	for (int i = mode + 1; i <= x; ++i) {
	    r *= rfunction(n1, n2, m1, psi, i);
	}
    }
    else if (x < mode) {
	for (int i = mode - 1; i >= x; --i) {
	    r /= rfunction(n1, n2, m1, psi, i + 1);
	}
    }

    return r;
}

/*
 * Calculates the normalized density of x
 */
static double density_normalized(int x, int n1, int n2, int m1, double psi)
{
    int ll = max(0, m1 - n2);
    int uu = min(n1, m1);

    if (x < ll || x > uu) return 0;

    // Density at mode has reference value 1
    int mode = modeCompute(n1, n2, m1, psi);
    double psum = 1, px = 1;

    // Calculate relative density above the mode
    if (mode < uu) {
	double r = 1;
	for (int i = mode + 1; i <= uu; ++i) {
	    r *= rfunction(n1, n2, m1, psi, i);
	    if (i == x) px = r;
	    psum += r;
	}
    }
    // Calculate relative density below the node
    if (mode > ll) {
	double r = 1;
	for (int i = mode - 1; i >= ll; --i) {
	    r /= rfunction(n1, n2, m1, psi, i + 1);
	    if (i == x) px = r;
	    psum += r;
	}
    }

    // Normalize
    return px/psum;
}

/*
 * Returns a vector of normalized probabilities for a hypergeometric
 * distribution. If the returned vector p is of length N, then p[0]
 * represents the probability of the lowest possible value ll = max(0,
 * m1 - n2) and p[N - 1] represents the probability of the largest
 * possible value uu = min(n1, m1).
 */
static vector<double> density_full(int n1, int n2, int m1, double psi)
{
    int ll = max(0, m1 - n2);
    int uu = min(n1, m1);
    unsigned int N = static_cast<unsigned int>(uu - ll + 1);
    vector<double> p(N);

    // Density at mode has reference value 1
    int mode = modeCompute(n1, n2, m1, psi);
    p[mode - ll] = 1;
    double sump = 1;

    // Calculate relative density above the mode
    if (mode < uu) {
	double r = 1;
	for (int i = mode + 1; i <= uu; ++i) {
	    r *= rfunction(n1, n2, m1, psi, i);
	    p[i - ll] = r;
	    sump += r;
	}
    }
    // Calculate relative density below the node
    if (mode > ll) {
	double r = 1;
	for (int i = mode - 1; i >= ll; --i) {
	    r /= rfunction(n1, n2, m1, psi, i + 1);
	    p[i - ll] = r;
	    sump += r;
	}
    }
  
    //Normalize
    for (int i = 0; i < N; ++i) {
	p[i] /= sump;
    }
    return p;
}

/*
 * Sample from a discrete unimodal distribution given a vector of
 * normalized probabilities pi[] and an estimate of the mode. U should
 * be drawn from a uniform(0,1) distribution. We spend our probability
 * allocation U as rapidly as possible, starting at the mode and then
 * stepping out to the left and the right.
 * 
 * This function is quite general and could be used for other discrete
 * distributions when we have a good guess of the mode. The function
 * also works if "mode" is not actually the mode, but is
 * correspondingly less efficient.
 */
static int sampleWithMode(int mode, vector<double> const &pi, double U) 
{
    int N = pi.size();
    int left = mode - 1;
    int right = mode;
    while (left >= 0 || right < N) {
	if (right < N && (left < 0 || pi[right] > pi[left])) {
	    U -= pi[right];
	    if (U <= 0) return right;
	    else ++right;
	}
	else {
	    U -= pi[left];
	    if (U <= 0) return left;
	    else --left;
	}
    }

    //Only reachable with bad arguments
    return mode;
}

double DHyper::d(double z, PDFType type,
                 vector<double const *> const &parameters, 
		 bool give_log) const
{
    int n1,n2,m1;
    double psi;
    getParameters(n1, n2, m1, psi, parameters);

    int x = static_cast<int>(z);
    double den = 0;
    if (type == PDF_PRIOR) {
	den = density_unnormalized(x, n1, n2, m1, psi);
    }
    else {
	den = density_normalized(x, n1, n2, m1, psi);
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

    int ll = max(0, m1 - n2);
    int uu = min(n1, m1);

    double sumpi = 0;
    if (x >= ll) {
	if (x >= uu) {
	    sumpi = 1;
	}
	else {
	    vector<double> pi = density_full(n1, n2, m1, psi);
	    for (int i = ll; i <= x; ++i) {
		sumpi += pi[i - ll];
	    }
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
    double tol = 64 * DBL_EPSILON;
    int n1,n2,m1;
    double psi;
    getParameters(n1, n2, m1, psi, parameters);
  
    int ll = max<int>(0, m1 - n2);
    int uu = min(n1, m1);
  
    vector<double> pi = density_full(n1, n2, m1, psi);

    if (log_p)
	p = exp(p);
    if (!lower)
	p = 1 - p;

    double sumpi = 0;
    for (int i = ll; i < uu; ++i) {
	sumpi += pi[i - ll];
	if (sumpi > p - tol) {
	    return i;
	}
    }
    return uu;
}

double DHyper::r(vector<double const *> const &parameters, RNG *rng) const
{
    int n1,n2,m1;
    double psi;
    getParameters(n1, n2, m1, psi, parameters);

    int ll = max(0, m1 - n2);
    int mode = modeCompute(n1, n2, m1, psi);
    vector<double> pi = density_full(n1, n2, m1, psi);
    return sampleWithMode(mode - ll, pi, rng->uniform()) + ll;
}

double DHyper::l(vector<double const *> const &parameters) const
{
    int n1,n2,m1;
    double psi;
    getParameters(n1, n2, m1, psi, parameters);
  
    return max(0, m1 - n2);
}

double DHyper::u(vector<double const *> const &parameters) const
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

double DHyper::KL(vector<double const *> const &para,
		  vector<double const *> const &parb) const
{
    int n1a,n2a,m1a;
    double psia;
    getParameters(n1a, n2a, m1a, psia, para);

    int lla = max(0, m1a - n2a);
    int uua = min(n1a, m1a);

    int n1b,n2b,m1b;
    double psib;
    getParameters(n1b, n2b, m1b, psib, parb);

    int llb = max(0, m1b - n2b);
    int uub = min(n1b, m1b);

    if (lla < llb || uua > uub)
	return JAGS_POSINF;

    vector<double> da = density_full(n1a, n2a, m1a, psia);
    vector<double> db = density_full(n1b, n2b, m1b, psib);

    double y = 0;
    for (int i = lla; i <= uua; ++i) {
	double proba = da[i - lla];
	double probb = db[i - llb];
	y += proba * (log(proba) - log(probb));
    }
    return y;
}

}}
