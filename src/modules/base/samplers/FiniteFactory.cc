#include <config.h>
#include <sampler/ImmutableSampler.h>
#include <sampler/SingletonGraphView.h>

#include "FiniteMethod.h"
#include "FiniteFactory.h"

#include <vector>

using std::vector;
using std::string;

namespace jags {
namespace base {

    bool 
    FiniteFactory::canSample(StochasticNode *snode, Graph const &) const
    {
	return FiniteMethod::canSample(snode);
    }

    Sampler * FiniteFactory::makeSampler(StochasticNode *snode,
					 Graph const &graph) const
    {
	SingletonGraphView *gv = new SingletonGraphView(snode, graph);
	FiniteMethod *method = new FiniteMethod(gv);
	return new ImmutableSampler(gv, method, name());
    }

    string FiniteFactory::name() const
    {
	return "base::Finite";
    }

}}
