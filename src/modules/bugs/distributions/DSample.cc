#include <config.h>
#include "DSample.h"
#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <module/ModuleError.h>

#include <cmath>
#include <algorithm>
#include <list>
#include <numeric>

using std::vector;
using std::list;
using std::log;
using std::min;
using std::fill;
using std::accumulate;
using std::sort;
using std::exp;

static inline double  SIZE(vector<double const *> const &par) {
    return *par[1];
}

/* Utility functions used to sort probabilities in descending order */

static inline bool gt_doubleptr (double const *arg1, double const *arg2) {
  return *arg1 > *arg2;
}

static inline bool gt_double (double arg1, double arg2) {
  return arg1 > arg2;
}

namespace jags {
    namespace bugs {

	DSample::DSample()
	    : VectorDist("dsample", 2) 
	{}
	
	bool DSample::isDiscreteValued(vector<bool> const &) const
	{
	    return true; //Takes only 0, 1 values
	}

	bool
	DSample::checkParameterLength(vector<unsigned long> const &len) const
	{
	    //Check that PROB is non-empty and SIZE is a scalar
	    return len[0] >= 1 && len[1] == 1;
	}

	bool DSample::checkParameterDiscrete(vector<bool> const &mask) const
	{
	    return mask[1]; //SIZE is discrete-valued
	}

	bool DSample::checkParameterValue(vector<double const *> const &par,
					vector<unsigned long> const &len) const
	{
	    if (SIZE(par) < 0 || SIZE(par) > len[0]) {
		//Too few or too many samples requested
		return false;
	    }
	    
	    for (unsigned long i = 0; i < len[0]; ++i) {
		if (par[0][i] <= 0) {
		    //Zero probablity of being sampled
		    return false;
		}
	    }
	    return true;
	}

	unsigned long DSample::length(vector<unsigned long> const &len) const
	{
	    //Length of sample is the same as length of probability weights
	    return len[0];
	}

	double
	DSample::logDensity(double const *x, PDFType type,
			    vector<double const *> const &par,
			    vector<unsigned long> const &parlen,
			    double const *lower, double const *upper) const
	{
	    /* Basic sanity checks: is x consistent with parameters */

	    unsigned long T = parlen[0]; // Length of x
	    unsigned long K = 0; // Count number of sampled values
	    for (unsigned long t = 0; t < T; ++t) {
		if (x[t] == 1) {
		    ++K;
		}
		else if (x[t] != 0) {
		    return JAGS_NEGINF;
		}
	    }
	    if (SIZE(par) != K) {
		return JAGS_NEGINF;
	    }

	    /* If there are more sampled than unsampled values then
	     * we turn the calculations around: counting unsampled values
	     * and inverting the probability weights.
	     */
	    int sign = 1;
	    int y = 1;
	    if (K > T/2) {
		sign = -1;
		y = 0;
		K = T - K;
	    }

	    /* Probability weights are normalized to have maximum zero on the
	     * log scale */
	    double const *p = par[0]; // Vector of probability weights
	    vector<double> logp(T); // Log probability weights

	    double lpmax = JAGS_NEGINF;
	    for (unsigned long t = 0; t < T; ++t) {
		logp[t] = sign * log(p[t]);
		if (logp[t] > lpmax) {
		    lpmax = logp[t];
		}
	    }

	    /* Likelihood calculations */
	    
	    double loglik = 0;

	    /* Numerator: contribution to log likelihood of all
	     * sampled values */
	    for (unsigned long t = 0; t < T; ++t) {
		if (x[t] == y) {
		    loglik += logp[t] - lpmax;
		}
	    }

	    if (type != PDF_PRIOR) {
	    
		/* Denominator: Recursively calculate contributions over
		 * all possible case sets of size K */
		
		/* 
		   For numerical stability, sort elements of logp in
		   descending order. This ensures that we are adding
		   increasingly smaller terms to f[k] for k in 1:(K+1)
		*/
		sort(logp.begin(), logp.end(), gt_double);
		vector<double> f(K + 1, 0);
		f[0] = 1;
		for (unsigned long t = 0; t < T; ++t) {
		    double Ct = exp(logp[t] - lpmax);
		    for (unsigned long k = min<unsigned long>(K, t+1); k > 0; --k)
		    {
			f[k] += Ct * f[k-1];
		    }
		}
		loglik -= log(f[K]);
	    }

	    return loglik;
	}


	void DSample::randomSample(double *x,
				   vector<double const *> const &par,
				   vector<unsigned long> const &parlen,
				   double const *, double const *,
				   RNG *rng) const
	{
	    double const * const probs = par[0];
	    const unsigned long N = parlen[0];
	    
	    //Create a vector of pointers to the elements of the vector
	    //of probability weights. Sort them in reverse order.
	    list<double const *> pptrs(N);
	    list<double const *>::iterator q;
	    double const *y;
	    for (q = pptrs.begin(), y = par[0]; q != pptrs.end(); ++q) {
		*q = y++;
	    }
	    pptrs.sort(gt_doubleptr);

	    //Initialize sample value by setting all elements to zero
	    fill(x, x + N, 0);

	    double sump = accumulate(probs, probs + N, 0.0);

	    unsigned long K = static_cast<unsigned long>(SIZE(par));
	    for (unsigned long k = 0; k < K; ++k) {
		double prand = sump * rng->uniform();
		for (q = pptrs.begin(); q != pptrs.end(); ++q) {
		    prand -= **q;
		    if (prand <= 0) {
			unsigned long i = *q - probs;
			x[i] = 1;
			sump -= **q;
			pptrs.erase(q);
			break;
		    }
		}

	    }
	}

	void DSample::support(double *lower, double *upper,
			    vector<double const *> const &par,
			    vector<unsigned long> const &parlen) const
	{
	    for (unsigned long i = 0; i < parlen[0]; ++i) {
		lower[i] = 0;
		upper[i] = 1;
	    }
	}

	bool DSample::isSupportFixed(vector<bool> const &) const
	{
	    return true;
	}
	
	unsigned long DSample::df(vector<unsigned long> const &parlen) const
	{
	    return parlen[0];
	} 

    }
}
