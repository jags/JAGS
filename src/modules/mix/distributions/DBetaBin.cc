#include <config.h>
#include "DBetaBin.h"

#include <algorithm>

#include <JRmath.h>
#include <cfloat>

using std::vector;
using std::max;
using std::min;
using std::string;

#define A(par) (*par[0])
#define B(par) (*par[1])
#define SIZE(par) (*par[2])

static inline double dbb(double x, double a, double b, double n)
{
    //Probability density function
    //return choose(a+x-1, x) * choose(b+n-x-1, n-x) / choose(a+b+n-1, n);
    return choose(n, x) * beta(x + a, n - x + b) / beta(a, b);
}

static inline double ldbb(double x, double a, double b, double n)
{
    //Log probability density function
    return lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b);
}

namespace jags {
namespace mix {

DBetaBin::DBetaBin()
    : RScalarDist("dbetabin", 3, DIST_SPECIAL, true)
{}

string DBetaBin::alias() const
{
    return "dbetabinom";
}

bool 
DBetaBin::checkParameterDiscrete (vector<bool> const &mask) const
{
    return mask[2] == true;
}

bool DBetaBin::checkParameterValue (vector<double const *> const &par) const
{
    return (SIZE(par) >= 0 && A(par) > 0.0 && B(par) > 0.0);
}

double DBetaBin::d(double x, PDFType type, vector<double const *> const &par, 
		   bool give_log) const
{
    if (give_log)
	return ldbb(x, A(par), B(par), SIZE(par));
    else
	return dbb(x, A(par), B(par), SIZE(par));
}
    
static double pbb(double x, double a, double b, double n)
{
    //Distribution function

    if (x < 0) return 0;
    if (x >= n) return 1;
    
    double y = 0;
    for (int i = 0; i <= x; ++i) {
	y += dbb(i, a, b, n);
    }
    return y;
}

double DBetaBin::p(double x, vector<double const *> const &par, 
		   bool lower, bool give_log) const
{
    double y = pbb(x, A(par), B(par), SIZE(par));
    if (!lower) y = 1 - y;
    return give_log ? log(y) : y;
}

static double qbb(double p, double a, double b, double n)
{
    //Quantile function
    
    if (p < 0) return 0;
    if (p >= 1) return n;

    /* fuzz to ensure left continuity: */
    p *= 1 - 64*DBL_EPSILON;
	    
    double psum = 0;
    for (int i = 0; i < n; ++i) {
	psum += dbb(i, a, b, n);
	if (psum > p) return i;
    }
    return n;
}

double DBetaBin::q(double p, vector<double const *> const &par, 
		   bool lower, bool log_p) const
{
    if (log_p) p = exp(p);
    if (!lower) p = 1 - p;
    return qbb(p, A(par), B(par), SIZE(par));
}

double DBetaBin::r(vector<double const *> const &par, RNG *rng) const
{
    return rbinom(SIZE(par), rbeta(A(par), B(par), rng), rng);
}

double DBetaBin::l(vector<double const *> const &par) const
{
    return 0;
}

double DBetaBin::u(vector<double const *> const &par) const
{
    return SIZE(par);
}

bool DBetaBin::isSupportFixed(vector<bool> const &fixmask) const
{
   return fixmask[2]; //SIZE is fixed;
} 

}}
