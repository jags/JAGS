#include "LMDesign.h"

#include <sampler/SingletonGraphView.h>
#include <sampler/Linear.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>

using std::vector;

namespace jags {    
    namespace bugs {

	static vector<Node const *> mkLinearPredictor(vector<StochasticNode *> const &stoch_children)
	{
	    vector<Node const *> lpnodes;
	    for (auto p = stoch_children.begin(); p != stoch_children.end(); ++p) {
		Node const *lp = (*p)->parents()[0];
		lpnodes.push_back(lp);
	    }
	    return lpnodes;
	}
	
	LMDesign::LMDesign(SingletonGraphView const *gv)
	    : SingletonForwardDiff(gv, mkLinearPredictor(gv->stochasticChildren()))
	{
	    if (!checkLinear(gv, false, false)) {
		throwLogicError("Error in LMDesign: GraphView is not a linear model");
	    }

	    vector<StochasticNode *> const &stoch_children = gv->stochasticChildren();
	    for (auto p = stoch_children.begin(); p != stoch_children.end(); ++p) {
		vector<Node const *> const &pparents = (*p)->parents();
		for (unsigned long i = 1; i < pparents.size(); ++i) {
		    if (gv->isDependent(pparents[i])) {
			throwLogicError("Error in LMDesign: GraphView is not a linear model");
		    }
		}
	    }
	}
	
    } /* Namespace bugs */
} /* Namespace jags */
