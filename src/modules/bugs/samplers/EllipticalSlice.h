/* Copyright (C) 2026 Marcel Jonker */

#ifndef ELLIPTICAL_SLICE_H_
#define ELLIPTICAL_SLICE_H_

#include <sampler/MutableSampleMethod.h>
#include <vector>

namespace jags {

class SingletonGraphView;

    namespace bugs {

	/**
	 * Elliptical slice sampler for an unbounded multivariate normal
	 * stochastic node with a non-conjugate likelihood.
	 *
	 * The sampler follows Murray, Adams and MacKay (2010), "Elliptical
	 * slice sampling", AISTATS. It has no tunable parameters and does not
	 * require an adaptation phase.
	 */
	class EllipticalSlice : public MutableSampleMethod
	{
	    SingletonGraphView const *_gv;
	    unsigned int _chain;
	    unsigned long _length;
	    std::vector<double> _xcur;
	    std::vector<double> _xprop;
	    std::vector<double> _nu;
	    std::vector<double> _chol;
	    bool _fixed_prec;

	public:
	    /**
	     * @param gv GraphView object wrapping a single dmnorm node
	     *
	     * @param chain Chain number, starting from zero
	     */
	    EllipticalSlice(SingletonGraphView const *gv, unsigned int chain);
	    void update(RNG *rng) override;
	    bool isAdaptive() const override;
	    void adaptOff() override;
	    bool checkAdaptation() const override;
	};

    }
}

#endif /* ELLIPTICAL_SLICE_H_ */
