#include <config.h>

#include "ScaledWishartFactory.h"
#include "ScaledWishart.h"

#include <sampler/MutableSampler.h>
#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>

using std::vector;
using std::string;

namespace jags {
    namespace glm {

	ScaledWishartFactory::~ScaledWishartFactory()
	{
	    //Nothing to do here; only for vtable.
	}
	
	bool ScaledWishartFactory::canSample(StochasticNode *snode,
					   Graph const &graph) const
	{
	    return ScaledWishart::canSample(snode, graph);
	}
  
	Sampler *ScaledWishartFactory::makeSampler(StochasticNode *snode, 
						Graph const &graph) const
	{
	    unsigned int nchain = snode->nchain();
	    vector<MutableSampleMethod*> methods(nchain, 0);
    
	    SingletonGraphView *gv = new SingletonGraphView(snode, graph);
	
	    for (unsigned int ch = 0; ch < nchain; ++ch) {
		methods[ch] = new ScaledWishart(gv, ch);
	    }
    
	    return new MutableSampler(gv, methods, "glm::ScaledWishart");
	}

	string ScaledWishartFactory::name() const
	{
	    return "glm::ScaledWishart";
	}
  
    }
}
