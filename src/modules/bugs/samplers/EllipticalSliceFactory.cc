/* Copyright (C) 2026 Marcel Jonker */

#include <config.h>

#include "EllipticalSlice.h"

#include "EllipticalSliceFactory.h"
#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>
#include <sampler/MutableSampler.h>
#include <sampler/SingletonGraphView.h>
#include <string>
#include <vector>

using std::vector;
using std::string;

namespace jags {
    namespace bugs {

	bool
	EllipticalSliceFactory::canSample(StochasticNode * snode, Graph const &) const
	{
	    /* The sampler sets the value of the whole node, so it cannot be
	       used for partially observed nodes. */
	    return snode->distribution()->name() == "dmnorm" && !isBounded(snode)
		&& !isObserved(snode);
	}

	Sampler *
	EllipticalSliceFactory::makeSampler(StochasticNode *snode, Graph const &graph) const
	{
	    unsigned int N = snode->nchain();
	    vector<MutableSampleMethod*> methods(N, nullptr);

	    SingletonGraphView *gv = new SingletonGraphView(snode, graph);
	    for (unsigned int ch = 0; ch < N; ++ch) {
		methods[ch] = new EllipticalSlice(gv, ch);
	    }
	    return new MutableSampler(gv, methods, "bugs::EllipticalSlice");
	}

	string
	EllipticalSliceFactory::name() const
	{
	    return "bugs::EllipticalSlice";
	}

    }
}
