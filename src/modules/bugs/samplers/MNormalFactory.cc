#include <config.h>

#include "MNormal.h"
#include "MNormalFactory.h"
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
MNormalFactory::canSample(StochasticNode * snode, Graph const &) const
{
    return snode->distribution()->name() == "dmnorm" && !isBounded(snode);
}

Sampler *
MNormalFactory::makeSampler(StochasticNode *snode, Graph const &graph) const
{
    unsigned int N = snode->nchain();
    vector<MutableSampleMethod*> methods(N, nullptr);

    SingletonGraphView *gv = new SingletonGraphView(snode, graph);
    for (unsigned int ch = 0; ch < N; ++ch) {
        methods[ch] = new MNormMetropolis(gv, ch);
    }
    return new MutableSampler(gv, methods, "bugs::MNormalMetropolis");
}

string MNormalFactory::name() const 
{
    return "bugs::MNormal";
}

}}
