#include <config.h>
#include <distribution/Distribution.h>
#include <graph/Graph.h>
#include <graph/StochasticNode.h>
#include <graph/NodeError.h>
#include <sampler/MutableSampler.h>
#include <sampler/GraphView.h>

#include "DSumFactory.h"
#include "RealDSum.h"
#include "DiscreteDSum.h"
#include "DMultiDSum.h"

#include <algorithm>

using std::list;
using std::vector;
using std::string;

namespace jags {
namespace bugs {

static StochasticNode const *getDSumChild(StochasticNode *node)
{
    list<StochasticNode*>::const_iterator p;
    for (p = node->stochasticChildren()->begin(); 
	 p != node->stochasticChildren()->end(); ++p) 
    {
	//Skip unobserved nodes
	if (isObserved(*p) && (*p)->distribution()->name() == "dsum") {
	    if (isParameter(*p)) {
		//Node is only partly observed. We can't handle this case
		return nullptr;
	    }
	    else {
		return *p;
	    }
	}
	    
    }
    return nullptr;
}

Sampler * DSumFactory::makeSampler(list<StochasticNode*> const &nodes,
				   Graph const &graph) const
{
    //Find DSum node
    StochasticNode const *dsum_node = nullptr;
    for (list<StochasticNode*>::const_iterator p = nodes.begin(); 
	 p != nodes.end(); ++p) 
    {
	dsum_node = getDSumChild(*p);
	if (dsum_node)
	    break;
    }
    
    if (!dsum_node)
	return nullptr;
    
    //See if we can sample the parents. This can only be done if they
    //are unobserved stochastic nodes in the sample set

    vector<StochasticNode *> parameters;
    vector<Node const *> const &parents = dsum_node->parents();
    vector<Node const *>::const_iterator pp;
    for (pp = parents.begin(); pp != parents.end(); ++pp) {
	list<StochasticNode*>::const_iterator q =
	    find(nodes.begin(), nodes.end(), *pp);
	if (q != nodes.end()) {
	    parameters.push_back(*q);
	}
	else {
	    return nullptr;
	}
    }
    
    bool discrete;
    bool multinom = false;
    string name;
    if (RWDSum::canSample(parameters, graph, false, false)) {
	discrete = false;
	name = "bugs::RealDSum";
    }
    else if (RWDSum::canSample(parameters, graph, true, false)) {
	discrete = true;
	name = "bugs::DiscreteDSum";
    }
    else if (RWDSum::canSample(parameters, graph, true, true)) {
	discrete = true;
	multinom = true;
	name = "bugs::DMultiDSum";
    }
    else {
	return nullptr;
    }
    
    GraphView *gv = new GraphView(parameters, graph, true);
    unsigned int nchain = parameters[0]->nchain();
    vector<MutableSampleMethod*> methods(nchain, nullptr);
    for (unsigned int ch = 0; ch < nchain; ++ch) {
	if (discrete) {
	    if (multinom)
		methods[ch] = new DMultiDSum(gv, ch);
	    else
		methods[ch] = new DiscreteDSum(gv, ch);
	}
	else {
	    methods[ch] = new RealDSum(gv, ch);
	}
    }
    return new MutableSampler(gv, methods, name);

}

string DSumFactory::name() const
{
    return "bugs::DSum";
}

vector<Sampler*>  DSumFactory::makeSamplers(list<StochasticNode*> const &nodes, 
					    Graph const &graph) const
{
    Sampler *s = makeSampler(nodes, graph);
    if (s)
	return vector<Sampler*>(1, s);
    else 
	return vector<Sampler*>();
}

}}
