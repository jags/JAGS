#include <config.h>
#include <sarray/SArray.h>
#include "DNormMix.h"

#include <cmath>

#include <Rmath.h>

using std::vector;

static inline double const * MU(vector<SArray const *> const &par)
{
    return par[0]->value();
}

static inline double const * TAU(vector<SArray const *> const &par)
{
    return par[1]->value();
}

static inline double const * PROB(vector<SArray const *> const &par)
{
    return par[2]->value();
}

static inline long NCAT(vector<SArray const *> const &par)
{
    return par[0]->length();
}

DNormMix::DNormMix()
  : DistReal("dnormmix", 3, DIST_UNBOUNDED, false)
{}

DNormMix::~DNormMix()
{}

bool DNormMix::checkParameterDim(vector<Index> const &dims) const
{
    if (dims[0][0] == 1)
        return false; // Must be a mixture
   
    // Parameter dimensions must match (but they need not be vectors)
    return (dims[1] == dims[0]) && (dims[2] == dims[0]);
}

bool DNormMix::checkParameterValue(vector<SArray const *> const &par) const
{
    double const *tau = TAU(par);
    double const *prob = PROB(par);
    unsigned long Ncat = NCAT(par);
    double sump = 0.0;
    for (unsigned int i = 0; i < Ncat; ++i) {
	if (tau[i] <= 0)
	    return false;
	sump += prob[i];
    }
    return (fabs(sump - 1) < 16 * DBL_EPSILON);
}

double
DNormMix::d(double x, vector<SArray const *> const &par, bool give_log) const
{
    double const *mu = MU(par);
    double const *tau = TAU(par);
    double const *p = PROB(par);
    unsigned long Ncat = NCAT(par);

    double density = 0.0;
    for (unsigned int i = 0; i < Ncat; ++i) {
	density += p[i] * dnorm(x, mu[i], 1/sqrt(tau[i]), 0);
    }
    return give_log ? log(density) : density;
}

double
DNormMix::p(double q, vector<SArray const *> const &par, bool lower, bool give_log)
    const
{
    double const *mu = MU(par);
    double const *tau = TAU(par);
    double const *prob = PROB(par);
    unsigned long Ncat = NCAT(par);
    
    double density = 0.0;
    for (unsigned int i = 0; i < Ncat; ++i) {
	density += prob[i] * pnorm(q, mu[i], 1/sqrt(tau[i]), lower, 0);
    }
    return give_log ? log(density) : density;
}

static double 
pmix(double q, double const *mu, double const *sigma, double const *p,
     unsigned int Ncat) 
{
    /* Probability distribution function of a normal mixture */

    double density = 0.0;
    for (unsigned int i = 0; i < Ncat; ++i) {
	density += p[i] * pnorm(q, mu[i], sigma[i], true, false);
    }
    return density;
}

static double 
qmix (double p, double const *mu, double const *sigma,
      double const *prob, unsigned int Ncat)
{
    /* Calculate quantile of normal mixture by bisection
       Based on an algorithm in the R package ensembleBMA by
       Adrian E. Raftery, J. McLean Sloughter, and Michael Polakowski, 
    */

    const double C = qnorm(p, 0, 1, true, false);
    double lower = mu[0] - C * sigma[0];
    double upper = mu[0] + C * sigma[0];
    double mid = (lower + upper)/2;

    while (upper - lower > 16 * sqrt(DBL_EPSILON)) {
	if (pmix(mid, mu, sigma, prob, Ncat) > p) {
	    upper = mid;
        }
        else {
            lower = mid;
        }
        mid = (lower + upper)/2;
    }
    return mid;
}

double 
DNormMix::q(double p, vector<SArray const *> const &par, bool lower, bool log_p)
    const
{
    if (!lower)
	p = 1 - p;
    if (log_p)
	p = exp(p);

    double const *mu = MU(par);
    double const *tau = TAU(par);
    double const *prob = PROB(par);
    unsigned long Ncat = par[0]->length();
    
    double * sigma = new double[Ncat];
    for (unsigned int i = 0; i < Ncat; ++i) {
	sigma[i] = 1/sqrt(tau[i]);
    }
    double q = qmix(p, mu, sigma, prob, Ncat);
    delete [] sigma;

    return q;
}

double 
DNormMix::r(vector<SArray const *> const &par, RNG *rng) const
{
    double const *mu = MU(par);
    double const *tau = TAU(par);
    double const *p = PROB(par);
    unsigned long Ncat = par[0]->length();    
    
    // Select mixture component (r)
    unsigned int r = Ncat - 1;
    double sump = 0;
    double p_rand = runif(0, 1, rng);
    for (unsigned int i = 0; i < Ncat - 1; ++i) {
	sump += p[i];
	if (sump > p_rand) {
	    r = i;
	    break;
	}
    }

    // Now sample from conditional distribution of component r
    double ans = rnorm(mu[r], 1/sqrt(tau[r]), rng);

    return ans;
}

void 
DNormMix::typicalValue(SArray& x, std::vector<SArray const *> const &par) const
{
    double const *mu = MU(par);
    double const *prob = PROB(par);
    unsigned long Ncat = NCAT(par);

    unsigned long j = 0;
    for (unsigned long i = 1; i < Ncat; ++i) {
       if (prob[i] > prob[j]) {
           j = i;
       }
    }
    x.setValue(mu[j], 0);
}
