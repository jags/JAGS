#include "WAICMonitorFactory.h"
#include "WAICMonitor.h"

#include <model/BUGSModel.h>
#include <graph/StochasticNode.h>

using std::string;
using std::vector;

namespace jags {
    namespace dic {
    
	Monitor *WAICMonitorFactory::getMonitor(string const &name,
						Range const &range,
						BUGSModel *model,
						string const &type,
						string &msg)
	{
	    if (type != "mean")
		return 0;
	
	    if (name != "WAIC")
		return 0;

	    if (!isNULL(range))
		msg = string("Cannot monitor a subset of ") + name;
	
	    vector<StochasticNode const *> observed_nodes;
	    vector<StochasticNode *> const &snodes = model->stochasticNodes();
	    for (vector<StochasticNode *>::const_iterator p = snodes.begin();
		 p != snodes.end(); ++p)
	    {
		if (isObserved(*p)) {
		    observed_nodes.push_back(*p);
		}
	    }
	    if (observed_nodes.empty()) {
		msg = "There are no observed stochastic nodes";
		return 0;
	    }

	    Monitor *m = new WAICMonitor(observed_nodes);
	    m->setName(name);
	    vector<string> onames(observed_nodes.size());
	    for (unsigned int i = 0; i < observed_nodes.size(); ++i) {
		onames[i] = model->symtab().getName(observed_nodes[i]);
	    }
	    m->setElementNames(onames);

	    return m;
	}

	string WAICMonitorFactory::name() const
	{
	    return "dic::WAIC";
	}

    }
}
