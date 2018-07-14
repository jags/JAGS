#include <config.h>

#include "RW1.h"
#include "RW1Factory.h"

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
	RW1Factory::canSample(StochasticNode * snode, Graph const &) const
	{
	    return snode->distribution()->name() == "drw1";
	}

	Sampler *
	RW1Factory::makeSampler(StochasticNode *snode, Graph const &graph) const
	{
	    unsigned int N = snode->nchain();
	    vector<MutableSampleMethod*> methods(N, nullptr);
	    
	    SingletonGraphView *gv = new SingletonGraphView(snode, graph);
	    for (unsigned int ch = 0; ch < N; ++ch) {
		methods[ch] = new RW1(gv, ch);
	    }
	    return new MutableSampler(gv, methods, "bugs::RW1");
	}
	
	string RW1Factory::name() const 
	{
	    return "bugs::RW1";
	}
	
    }
}
