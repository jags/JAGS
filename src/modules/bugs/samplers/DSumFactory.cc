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
    set<StochasticNode*> dsum_nodes;
    set<StochasticNode*>::const_iterator p;

    //Find DSum nodes
    for (p = nodes.begin(); p != nodes.end(); ++p) {
	if ((*p)->distribution()->name() == "dsum")
	    dsum_nodes.insert(*p);
    }
  
    //See if we can sample them
    for (p = dsum_nodes.begin(); p != dsum_nodes.end(); ++p) {

	bool cansample = true;

	/* Get parents of dsumnode as a vector of Stochastic Nodes */
	vector<StochasticNode *> parameters;
	vector<Node const *> const &parents = (*p)->parents();
	vector<Node const *>::const_iterator pp;
	for (pp = parents.begin(); pp != parents.end(); ++pp) {
	    set<StochasticNode *>::const_iterator q =
		find(nodes.begin(), nodes.end(), *pp);
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
