#ifndef FINITE_METHOD_H_
#define FINITE_METHOD_H_

#include <sampler/ImmutableSampleMethod.h>

namespace jags {
    namespace base {

	/**
	 * Sampler for discrete distributions with support on a finite set.
	 */
	class FiniteMethod : public ImmutableSampleMethod {
	    SingletonGraphView const * const _gv;
	    const double _lower, _upper;
	    unsigned long _size;
	  public:
	    FiniteMethod(SingletonGraphView const *gv);
	    void update(unsigned int chain, RNG *rng) const;
	    static bool canSample(StochasticNode const *snode);
	};
	
    }
}

#endif /* FINITE_METHOD_H_ */
