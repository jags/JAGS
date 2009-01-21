#include <config.h>
#include <distribution/Distribution.h>
//asStochastic
#include <graph/StochasticNode.h>
#include <graph/Graph.h>
#include <graph/NodeError.h>
#include <sampler/DensitySampler.h>

#include "DSumFactory.h"
#include "DSumSampler.h"

#include <stdexcept>
#include <algorithm>

using std::set;
using std::vector;
using std::runtime_error;

void DSumFactory::makeSampler(set<StochasticNode*> &nodes,
			      Graph const &graph,
			      vector<Sampler*> &samplers) const
{
    set<Node*>::const_iterator p;
    vector<StochasticNode const*> dsum_nodes;

    //Find DSum nodes
    for (p = graph.nodes().begin(); p != graph.nodes().end(); ++p) {
	StochasticNode const *snode = asStochastic(*p);
	if (snode && snode->distribution()->name() == "dsum")
	    dsum_nodes.push_back(snode);
    }
  
    //See if we can sample them
    for (unsigned int i = 0; i < dsum_nodes.size(); ++i) {

	bool cansample = true;

	/* Get parents of dsumnode as a vector of Stochastic Nodes */
	vector<StochasticNode *> parameters;
	vector<Node const *> const &parents = dsum_nodes[i]->parents();
	for (unsigned int j = 0; j < parents.size(); ++j) {
	    set<StochasticNode *>::const_iterator q =
		find(nodes.begin(), nodes.end(), parents[j]);
	    if (q == nodes.end()) {
		cansample = false;
		break;
	    }
	    else {
		parameters.push_back(*q);
	    }
	}

	if (cansample && DSumMethod::canSample(parameters, graph)) {
	    for (unsigned int i = 0; i < parameters.size(); ++i) {
		nodes.erase(parameters[i]);
	    }
	    
	    unsigned int nchain = parameters[0]->nchain();
	    vector<DensityMethod*> methods(nchain, 0);
	    for (unsigned int ch = 0; ch < nchain; ++ch) {
		methods[ch] = new DSumMethod;
	    }
	    Sampler *sampler = new DensitySampler(parameters, graph, methods);
	    samplers.push_back(sampler);
	}
    }
}
