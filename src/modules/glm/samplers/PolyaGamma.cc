#include "PolyaGamma.h"
#include "Classify.h"
#include "PG.h"

#include <JRmath.h>
#include <graph/StochasticNode.h>
#include <rng/RNG.h>
#include <module/ModuleError.h>

#include <cmath>

using std::exp;
using std::log;
using std::sqrt;
using std::abs;

static const double one = 1;

namespace jags {
    namespace glm {

	//Optimal cutoff between left and right expansions of the
	//Jacobi density from Devroye
	static const double TRUNC=0.64;

	static inline double phi(double x) {
	    //Cumulative distribution function of a standard normal
	    return pnorm(x, 0.0, 1.0, true, false);
	}

	static double pigauss(double z)
	{
	    // P(X <= TRUNC) if X has an inverse Gaussian distribution with
	    // mu = 1/z; lambda = 1

	    const double ITRUNC_SQRT = 1 / sqrt(TRUNC);
	    double b = (TRUNC * z - 1) * ITRUNC_SQRT;
	    double a = -(TRUNC * z + 1) * ITRUNC_SQRT;
	    return phi(b) + exp(2 * z) * phi(a);
	}

	static double a(double n, double x)
	{
	    // Alternating series expansion of the Jacobi density
	    
	    n += 0.5;
	    double K = n * M_PI;
	    if (x > TRUNC) {
		return K * exp(-K * K * x / 2);
	    }
	    else if (x > 0) {
		double logy = -1.5 * log(x) - 2*n*n/x;
		return K * pow(M_2_PI, 1.5) * exp(logy);
	    }
	    else {
		return 0; //x == 0
	    }
	}


	static double rpolya_gamma(double z, RNG *rng)
	{
	    /* Sample from Polya-gamma PG(1, abs(z))
	    
	       In fact we are sampling from a Jacobi density, exploiting
	       the fact that PG(1, z) = J(1, z/2)/4; hence the transformation
	       of z on input and the return value on exit.
	    */
	    
	    z = abs(z)/2;

	    double K = M_PI * M_PI / 8 + z * z / 2;
	    double p = M_PI * exp(-K*TRUNC) / (2 * K);
	    double q = 2 * exp(-z) * pigauss(z);
	    
	    for (unsigned int i = 0; i < 10; ++i) {

		double X;
		if (rng->uniform() < p/(p+q)) {
		    // Sample from the tail with an exponential proposal
		    double E = rng->exponential();
		    X = TRUNC + E/K;
		}
		else {
		    // Sample from the body with a truncated IG proposal
		    X = rigauss(z, 1, TRUNC, rng);
		}
		double S = a(0, X);
		double Y = rng->uniform() * S;

		for(int n = 1; ; ++n) {
		    if (n % 2 == 1) {
			// odd terms
			S -= a(n, X);
			if (Y < S) return X/4;
		    }
		    else {
			// even terms
			S += a(n, X);
			if (Y > S) break;
		    }
		    if (n > 1000) {
			throwLogicError("Infinite loop in PolyaGamma?");
		    }
		}
	    }
	    throwLogicError("Failed to sample Polya-Gamma");
	    return 0; //-Wall
	} 

	static double const & getSize(StochasticNode const *snode,
				      unsigned int chain)
	{
	    if (getFamily(snode) == GLM_BERNOULLI) {
		return one;
	    }
	    else if (getFamily(snode) == GLM_BINOMIAL) {
		return *snode->parents()[1]->value(chain);
	    }
	    else {
		throwLogicError("Invalid outcome for PolyaGamma");
	    }
	    return one; //-Wall
	}
	
	PolyaGamma::PolyaGamma(StochasticNode const *snode, unsigned int chain)
	    : Outcome(snode, chain), _y(snode->value(chain)[0]),
	      _n(getSize(snode, chain)), _tau(1)
	{
	    //fixme: sanity checks on snode
	}
	
	double PolyaGamma::value() const 
	{
	    return (_y - _n/2)/_tau;
	}
	
	double PolyaGamma::precision() const 
	{
	    return _tau;
	}
	
	void PolyaGamma::update(RNG *rng)
	{
	    unsigned int N = static_cast<unsigned int>(_n);

	    _tau = 0.0;
	    for (unsigned int i = 0; i < N; ++i) {
		_tau += rpolya_gamma(_lp, rng);
	    }
	    if (_tau == 0.0) {
		throwLogicError("Bad Polya Gamma update");
	    }
	}

	bool PolyaGamma::canRepresent(StochasticNode const *snode)
	{
	    switch(getFamily(snode)) {
	    case GLM_BERNOULLI:
		break;
	    case GLM_BINOMIAL:
		if (!snode->parents()[1]->isFixed() ||
		    snode->parents()[1]->value(0)[0] > 19) {
		    /* 
		       PolyaGamma competes with AuxMixBinomial for
		       representing binomial outcomes with logit
		       link. PolyaGamma is more efficient for small N,
		       but its complexity is O(N), so it is much
		       slower for large N.

		       Where is the cutoff? For N = 20, the number of
		       mixture components used by AuxMixBinomial drops
		       from 9 to 4 (See LGMix.cc), so this may be a
		       good choice.
		    */
		    return false;
		}
		break;
	    default:
		return false;		
	    }
	    
	    return getLink(snode) == LNK_LOGIT;
	}
    }
}
