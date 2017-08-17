#include <config.h>

#include "ScaledGammaFactory.h"
#include "ScaledGamma.h"

#include <sampler/MutableSampler.h>
#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>

using std::vector;
using std::string;

namespace jags {
    namespace glm {

	ScaledGammaFactory::~ScaledGammaFactory()
	{
	    //Nothing to do here; only for vtable.
	}
	
	bool ScaledGammaFactory::canSample(StochasticNode *snode,
					   Graph const &graph) const
	{
	    return ScaledGamma::canSample(snode, graph);
	}
  
	Sampler *ScaledGammaFactory::makeSampler(StochasticNode *snode, 
						Graph const &graph) const
	{
	    unsigned int nchain = snode->nchain();
	    vector<MutableSampleMethod*> methods(nchain, 0);
    
	    SingletonGraphView *gv = new SingletonGraphView(snode, graph);
	
	    for (unsigned int ch = 0; ch < nchain; ++ch) {
		methods[ch] = new ScaledGamma(gv, ch);
	    }
    
	    return new MutableSampler(gv, methods, "glm::ScaledGamma");
	}

	string ScaledGammaFactory::name() const
	{
	    return "glm::ScaledGamma";
	}
  
    }
}
