#include "DOrderedLogit.h"

#include <JRmath.h>

namespace jags {
    namespace glm {
	
	DOrderedLogit::DOrderedLogit()
	    : DOrdered("dordered.logit")
	{}
	
	double DOrderedLogit::r(double mu, RNG *rng) const
	{
	    return rlogis(mu, 1.0, rng);
	}

	double DOrderedLogit::p(double x, double mu, bool lower, bool give_log)
	    const
	{
	    return plogis(x, mu, 1.0, lower, give_log);
	}

    }
}


	    
    
