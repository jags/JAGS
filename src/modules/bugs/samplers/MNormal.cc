#include <config.h>
#include "MNormal.h"

#include <DMNorm.h>

#include <graph/StochasticNode.h>
#include <sampler/SingletonGraphView.h>
#include <rng/RNG.h>
#include <matrix/lapack.h>
#include <matrix/blas.h>
#include <matrix/matrix.h>
#include <JRmath.h>

#include <cmath>
#include <algorithm>

using std::vector;
using std::exp;
using std::sqrt;
using std::min;
using std::max;
using std::pow;

namespace jags {
    namespace bugs {

	static const int one = 1;
	
	// Target effective sample size 
	static double ESS(unsigned long t, double b, double n0) {
	    if (t <= 1) {
		return n0 + t;
	    } else {
		return n0 + 1 + b*(t-1);
	    }
	}
	
	// Returns learning rate for the sample mean and variance for
	// a given sample size t.
	static double solve_lambda(unsigned long t, double b, double n0) {
	    double S0 = ESS(t-1, b, n0);
	    double S1 = ESS(t, b, n0);
	    double delta = S1 - S0;
	    
	    return (1 + sqrt(S0*(1-delta)/S1))/(1 + S0);
	}

	// Probability density function for a proposal moving from x0 to x1
	// (or vice versa)
	static double dstep(vector<double> const &x0, vector<double> const &x1,
			    vector<double> const &sigma_chol, double s)
	{
	    vector<double> y(x0);
	    for (unsigned long i = 0; i < y.size(); ++i) {
		y[i] -= x1[i];
	    }

	    int d = x0.size();
	    jags_dtrsv("L", "N", "N", &d, sigma_chol.data(), &d, y.data(), &one);

	    double ld = 0.0;
	    for (int i = 0; i < d; ++i) {
		ld -= y[i]*y[i];
	    }
	    ld /= 2*s*s;
	    return exp(ld);
	}

	// Randomly generate a step of size s. This must be added to the current
	// value to make a proposal
	static void rstep(vector<double> &x, vector<double> const &sigma_chol,
			  double s, RNG *rng)
	{
	    for (unsigned long i = 0; i < x.size(); ++i) {
		x[i] = s * rng->normal();
	    }
	    int d = x.size();
	    jags_dtrmv("L", "N", "N", &d, sigma_chol.data(), &d, x.data(), &one);
	}

	//Forward declaration of pforward is required because it is called by pbackward
	static double pforward(vector<vector<double>> const &xlist, vector<double> const &logdensity,
			       vector<double> const &sigma_chol, vector<double> const &ss,
			       unsigned long start, unsigned long end);
	
	static double pbackward(vector<vector<double>> const &xlist, vector<double> const &logdensity,
				vector<double> const &sigma_chol, vector<double> const & ss,
				unsigned long start, unsigned long end)
	{
	    if (end >= start) {
		return -1;
	    }
	    
	    unsigned long stage = start - end;
	    double alpha = exp(logdensity[end] - logdensity[start]);
	    for (unsigned long i = 1; i < stage; ++i) {
		double pf = pforward(xlist, logdensity, sigma_chol, ss, end, end+i);
		double pb = pbackward(xlist, logdensity, sigma_chol, ss, start, start-i);
		if (pb == 1) {
		    /* The reverse path is invalid because we should have accepted a previous
		       proposal with probability 1. */
		    return 0.0;
		}
		double N = (1 - pf) * dstep(xlist[end+i], xlist[end], sigma_chol, ss[i-1]);
		double D = (1 - pb) * dstep(xlist[start-i], xlist[start], sigma_chol, ss[i-1]);
		alpha *= N / D;
	    }
	    
	    return min(1.0, alpha);
	}

	static double pforward(vector<vector<double>> const &xlist, vector<double> const &logdensity,
			       vector<double> const &sigma_chol, vector<double> const &ss,
			       unsigned long start, unsigned long end)
	{
	    if (start >= end) return -1;

	    unsigned long stage = end - start;
	    double alpha = exp(logdensity[end] - logdensity[start]);
	    for (unsigned long i = 1; i < stage; ++i) {
		double pb = pbackward(xlist, logdensity, sigma_chol, ss, end, end-i);
		double pf = pforward(xlist, logdensity, sigma_chol, ss, start, start+i);
		if (pf == 1) {
		    /* The forward path is invalid because we should have accepted a previous
		       proposal with probability 1. NB This should never happen in practice. */
		    return 0.0;
		}
		double N = (1 - pb) * dstep(xlist[end-i], xlist[end], sigma_chol, ss[i-1]);
		double D = (1 - pf) * dstep(xlist[start+i], xlist[start], sigma_chol, ss[i-1]);
		alpha *= N / D;
	    }
	    return min(1.0, alpha);
	}

	/*
	  Calculates acceptance probability for a sequence of proposals under delayed rejection
	  by recursively calling pforward and pbackward
	*/
	static double delayed_rejection(vector<vector<double>> const &xlist,
					vector<double> const &logdensity,
					vector<double> const &sigma_chol,
					vector<double> const &ss,
					unsigned long ntry)
	{
	    return pforward(xlist, logdensity, sigma_chol, ss, 0, ntry);
	}

	
	static double cal_delta(double d, double a) {
	    /*
	      Optimal scaling for the Robbins-Munro algorithm to find
	      the target step size according to Garthwait, Fan, and
	      Sisson (2016) https://doi.org/10.1080/03610926.2014.936562
	    */
	    double A = - qnorm5(a/2, 0.0, 1.0, 1, 0);
	    return (1 - 1.0/d) * (sqrt(M_2PI) * exp(A*A/2)/(2*A)) + 1/(d*a*(1.0 - a));
	}

	/*
	static vector<double> initValue(SingletonGraphView const *gv, 
					unsigned int chain)
	{
	    unsigned long d = gv->length();
	    vector<double> ivalue(d);
	    gv->getValue(ivalue, chain);
	    return ivalue;
	}
	*/
	
	static vector<double> initSigma(SingletonGraphView const *gv,
					double prior_scale)
	{
	    /* Initializes Variance matrix to be a diagonal matrix
	     * with diagonal elements equal to prior_scale */
	    unsigned long d = gv->length();
	    vector<double> Sigma(d*d, 0);
	    for (unsigned long i = 0; i < d; ++i) {
		Sigma[i + d * i] = prior_scale;
	    }
	    return Sigma;
	}
	
	MNormMetropolis::MNormMetropolis(SingletonGraphView const *gv, 
					 unsigned int chain,
					 double prior_scale, unsigned int n0, double ess_fraction,
					 double target_p)
	    : _gv(gv), _chain(chain), _mu(gv->length(), 0), _Sigma(initSigma(gv, prior_scale)),
	    _Sigma_chol(initSigma(gv, sqrt(prior_scale))), _b(ess_fraction),
	    _n0(n0 + gv->length()), _ptarget(target_p), _t(0), _lstep(0), _lstep_bar(0), _pmean(0),
	    _delta(cal_delta(gv->length(), target_p)), _adapt(true)
	{
	    gv->checkFinite(chain); //Check validity of initial values
	}
	
	void MNormMetropolis::update(RNG *rng)
	{
	    const unsigned long ntry = 4;
	    
	    /* Set up bookkeeping for the delayed rejection algorithm */
	    vector<vector<double>> xlist(ntry+1);
	    vector<double> ldvec(ntry+1), ss(ntry+1);

	    unsigned long d = _gv->length();
	    vector<double> x(d);
	    _gv->getValue(x, _chain);
	    
	    // Store data for first step
	    xlist[0] = x;
	    ldvec[0] = _gv->logFullConditional(_chain);
	    ss[0] = 2.38 * exp(_lstep) / sqrt(d);

	    vector<double> y(d);
	    double alpha0 = 0; // Acceptance probability of first proposal
	    /* Delayed rejection loop */
	    for (unsigned long k = 0; k < ntry; ++k) {

		// New proposal
		rstep(y, _Sigma_chol, ss[k], rng);
		for (unsigned int i = 0; i < d; ++i) {
		    y[i] += x[i];
		}
		_gv->setValue(y, _chain);

		xlist[k+1] = y;
		ldvec[k+1] = _gv->logFullConditional(_chain);
		ss[k+1] = 0.5 * ss[k];
		
		double alpha = delayed_rejection(xlist, ldvec, _Sigma_chol, ss, k+1);
		if (k == 0) {
		    //Store this value for rescaling
		    alpha0 = alpha;
		}

		/* Acceptance step */
		if (rng->uniform() <= alpha) {
		    // Accept and break out of delayed rejection loop
		    break;
		} else {
		    // Reject and return to initial value
		    _gv->setValue(x, _chain);
		}
	    }
	    
	    if (_adapt) {
		rescale(alpha0);
	    }
	}
	
	void MNormMetropolis::adaptOff()
	{
	    _adapt = false;
	    _lstep = _lstep_bar;
	}
	
	void MNormMetropolis::rescale(double p)
	{
	    _t++;

	    unsigned long d = _gv->length();
	    vector<double> x(d);
	    _gv->getValue(x, _chain);
	    
	    
	    // Get learning rate for updating shape
	    double lambda = solve_lambda(_t, _b, _n0);
	    
	    vector<double> dx0(d), dx1(d);
	    for (unsigned long i = 0; i < d; ++i) {
		dx0[i] = x[i] - _mu[i];
		_mu[i] += lambda * dx0[i];
		dx1[i] = x[i] - _mu[i];
	    }
	    for (unsigned long i = 0; i < d; ++i) {
		_Sigma[i + d * i] += lambda * (dx0[i]*dx1[i] - _Sigma[i + d*i]);
		for (unsigned long j = 0; j < i; ++j) {
		    double V = (dx0[i]*dx1[j] + dx0[j]*dx1[i])/2;
		    _Sigma[i + d*j] += lambda * (V - _Sigma[i + d*j]);
		    _Sigma[j + d*i] = _Sigma[i + d*j];
		}
	    }

	    //Get Cholesky decomposition
	    cholesky(_Sigma_chol.data(), _Sigma.data(), d);
	   
	    _lstep += _delta * (p - _ptarget) * pow(_n0 + _t, -0.75);
	    _lstep = max(0.0, _lstep);
	    /*
	      The running mean of the log step sizes is the Ruppert-Polyak estimate
	      But we can't use it until we stop adapting
	    */
	    _lstep_bar += (_lstep - _lstep_bar)/_t;

	    //Monitor average step size with forgetting weights
	    _pmean += lambda * (p - _pmean);
	}
	
	bool MNormMetropolis::checkAdaptation() const
	{
	    if (_t < 2000) {
		return false;
	    }
	    if (_lstep >= 0.05) {
		return abs(_pmean - _ptarget) <= 0.05;
	    }
	    return true;
	}

	bool MNormMetropolis::isAdaptive() const
	{
	    return true;
	}
	  
    }
}
