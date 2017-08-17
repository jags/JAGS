#ifndef DORDERED_PROBIT_H_
#define DORDERED_PROBIT_H_

#include "DOrdered.h"

namespace jags {
    namespace glm {
	
	class DOrderedProbit : public DOrdered
	{
	public:
	    DOrderedProbit();
	    double r(double mu, RNG *rng) const;
	    double p(double x, double mu, bool lower, bool give_log) const;
	};

    }
}

#endif /* DORDERED_PROBIT_H_ */

	    
    
