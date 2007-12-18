#include "PDMonitorFactory.h"
#include "PDMonitor.h"

#include <model/Model.h>
#include <graph/StochasticNode.h>

#include <set>

using std::set;
using std::string;
using std::vector;

namespace dic {

    Monitor *PDMonitorFactory::getMonitor(Node const *node, 
					  Model *model,
					  unsigned int start,
					  unsigned int thin, 
					  string const &type)
    {
	if (type != "pD" || node->nchain() < 2)
	    return 0;

	unsigned int nchain = model->nchain();
	vector<RNG*> rngs;
	for (unsigned int i = 0; i < nchain; ++i) {
	    rngs.push_back(model->rng(i));
	}

	StochasticNode const *snode = asStochastic(node);
	if (snode)
	    return new PDMonitor(snode, start, thin, rngs, 10);
	else 
	    return 0;
    }

    vector<Node const*> 
    PDMonitorFactory::defaultNodes(Model *model, string const &type) const
    {
	vector<Node const*> dnodes;

        if (type == "pD" && model->nchain() >= 2) {
	    set<Node*> const &nodes = model->graph().nodes();
	    set<Node*>::const_iterator p = nodes.begin(); 
	    for ( ; p != nodes.end(); ++p) {
		if ((*p)->isObserved() && asStochastic(*p)) {
		    dnodes.push_back(*p);
		}
	    }
	}
	return dnodes;
    }
	
}