#include <config.h>

#include "BinomSlicer.h"
#include "BinomSliceFactory.h"

#include <sampler/MutableSampler.h>
#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>

#include <vector>

using std::vector;
using std::string;

namespace jags {
    namespace bugs {

	bool BinomSliceFactory::canSample(StochasticNode * node,
					  Graph const &graph) const
	{
	    return BinomSlicer::canSample(node, graph);
	}
	
	Sampler *BinomSliceFactory::makeSampler(StochasticNode *snode,
					   Graph const &graph) const
	{
	    unsigned int nchain = snode->nchain();
	    vector<MutableSampleMethod*> methods(nchain, nullptr);

	    SingletonGraphView *gv = new SingletonGraphView(snode, graph);

	    for (unsigned int ch = 0; ch < nchain; ++ch) {
		methods[ch] = new BinomSlicer(gv, ch);
	    }

	    return new MutableSampler(gv, methods, "bugs::BinomSlicer");
	}

	string BinomSliceFactory::name() const
	{
	    return "bugs::BinomSlice";
	}

    }
}
