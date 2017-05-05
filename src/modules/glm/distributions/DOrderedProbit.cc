#include "DOrderedProbit.h"

#include <JRmath.h>

namespace jags {
    namespace glm {
	
	DOrderedProbit::DOrderedProbit()
	    : DOrdered("dordered.probit")
	{}
	
	double DOrderedProbit::r(double mu, RNG *rng) const
	{
	    return rnorm(mu, 1.0, rng);
	}

	double
	DOrderedProbit::p(double x, double mu, bool lower, bool give_log) const
	{
	    return pnorm(x, mu, 1.0, lower, give_log);
	}

    }
}


	    
    
