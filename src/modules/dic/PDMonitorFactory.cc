#include "PDMonitorFactory.h"
#include "DefaultPDMonitor.h"
#include "KLPDMonitor.h"
#include "KLTab.h"

#include <model/Model.h>
#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>

#include <set>

using std::set;
using std::string;
using std::vector;

namespace dic {

    static KL const *findKL(string const &name) {
	static KLTab _kltab;
	return _kltab.find(name);
    }

    Monitor *PDMonitorFactory::getMonitor(Node const *node, 
					  Model *model,
					  unsigned int start,
					  unsigned int thin, 
					  string const &type)
    {
	if (type != "pD" || node->nchain() < 2)
	    return 0;
	
	StochasticNode const *snode = asStochastic(node);
	if (!snode)
	    return 0;
	
	if (isSupportFixed(snode)) {
	    
	    KL const *kl = findKL(snode->distribution()->name());
	    if (kl) {
		return new KLPDMonitor(snode, start, thin, kl);
	    }
	}
	
	unsigned int nchain = model->nchain();
	vector<RNG*> rngs;
	for (unsigned int i = 0; i < nchain; ++i) {
	    rngs.push_back(model->rng(i));
	}
	return new DefaultPDMonitor(snode, start, thin, rngs, 10);	
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
