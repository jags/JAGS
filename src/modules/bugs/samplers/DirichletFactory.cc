#include <config.h>

#include "Dirichlet.h"
#include "DirichletFactory.h"
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
DirichletFactory::canSample(StochasticNode * snode, Graph const &graph) const
{
    return snode->distribution()->name() == "ddirch";
}

Sampler *
DirichletFactory::makeSampler(StochasticNode *snode, Graph const &graph) const
{
    unsigned int nchain = snode->nchain();
    vector<MutableSampleMethod*> methods(nchain, nullptr);
    vector<StochasticNode*> nodes(1, snode);
    SingletonGraphView *gv = new SingletonGraphView(snode, graph);
    for (unsigned int ch = 0; ch < nchain; ++ch) {
        methods[ch] = new DirchMetropolis(gv, ch);
    }
    return new MutableSampler(gv, methods, "bugs::DirichletMetropolis");
}

string DirichletFactory::name() const
{
    return "bugs::Dirichlet";
}

}}
