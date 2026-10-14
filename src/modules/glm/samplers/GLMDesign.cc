#include "GLMDesign.h"

#include <sampler/SingletonGraphView.h>
#include <sampler/Linear.h>
#include <graph/StochasticNode.h>
#include <graph/LinkNode.h>
#include <module/ModuleError.h>

using std::vector;
using std::set;

namespace jags {    
    namespace glm {

	static vector<Node const *> mkLinearPredictor(vector<StochasticNode *> const &outcomes)
	{
	    vector<Node const *> linear_predictor(outcomes.size());
	    auto q = linear_predictor.begin();
	    for (auto p = outcomes.begin(); p != outcomes.end(); ++p) {
		Node const *eta = (*p)->parents()[0];
		if (dynamic_cast<LinkNode const *>(eta) != nullptr) {
		    eta = eta->parents()[0];
		}
		*q = eta;
		++q;
	    }
	    return linear_predictor;
	}
	
	GLMDesign::GLMDesign(SingletonGraphView const *gv,
			     vector<StochasticNode *> const &outcomes)
	    : SingletonForwardDiff(gv, mkLinearPredictor(outcomes))
	{
	    if (!checkLinear(gv, false, true)) {
		throwLogicError("GLMDesign: GraphView is not a GLM");
	    }

	    set<StochasticNode*> scset;
	    for (auto p = outcomes.begin(); p != outcomes.end(); ++p) {
		if (!scset.insert(*p).second) {
		    throwLogicError("GLMDesign: repeated outcome");
		}
	    }
	    
	    for (auto p = gv->stochasticChildren().begin(); p != gv->stochasticChildren().end(); ++p) {
		if (scset.count(*p) == 0) {
		    throwLogicError("GLMDesign: Outcomes not a superset of the stochastic children");
		}
		vector<Node const *> const &pparents = (*p)->parents();
		for (unsigned long i = 1; i < pparents.size(); ++i) {
		    if (gv->isDependent(pparents[i])) {
			throwLogicError("GLMDesign: variance also depends on sampled nodes"); 
		    }
		}
	    }
	}
	
    } /* Namespace glm */
} /* Namespace jags */
