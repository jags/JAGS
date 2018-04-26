#ifndef DISCRETE_SLICER_H_
#define DISCRETE_SLICER_H_

#include <sampler/Slicer.h>
#include <graph/Graph.h>

namespace jags {
    
    class StochasticNode;
    class SingletonGraphView;

namespace base {
    
    /**
     * @short Discrete slice sampler
     *
     * Slice sampler for discrete valued distributions.
     *
     * To sample a discrete valued random variable X, we create an
     * auxiliary variable Y which is real valued, and for which p_[Y](y) =
     * p_[X](floor(y)).  Hence Y has a piecewise-constant density function
     * on the real line.  We sample Y and then set X = floor(Y).
     */
    class DiscreteSlicer : public Slicer {
	SingletonGraphView const *_gv;
	const unsigned int _chain;
	double _x;
    public:
	/**
	 * Constructor for Slice Sampler
	 * @param width Initial width of slice
	 * @param ndoubles Maximum number of doublings of slice 
	 */
	DiscreteSlicer(SingletonGraphView const *gv, unsigned int chain,
                       double width=2, unsigned int ndoubles = 10);
	void setValue(double x);
	double value() const;
	void getLimits(double *lower, double *upper) const;
	void update(RNG*);
	static bool canSample(StochasticNode const *node);
	double logDensity() const;
    };

}}

#endif /* DISCRETE_SLICER_H_ */
