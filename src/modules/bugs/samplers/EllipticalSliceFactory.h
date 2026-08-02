/*Copyright (C) 2026 Marcel Jonker */

#ifndef ELLIPTICAL_SLICE_FACTORY_H_
#define ELLIPTICAL_SLICE_FACTORY_H_

#include "EllipticalSlice.h"
#include <sampler/SingletonFactory.h>

namespace jags {
    namespace bugs {

	/**
	 * @short Factory object for the elliptical slice sampler
	 */
	class EllipticalSliceFactory : public SingletonFactory
	{
	public:
	    bool canSample(StochasticNode *snode, Graph const &graph) const override;
	    Sampler *makeSampler(StochasticNode *snode, Graph const &graph) const override;
	    std::string name() const override;
	};

    }
}

#endif /* ELLIPTICAL_SLICE_FACTORY_H_ */
