#include <config.h>

#include "REFactory2.h"
#include "REMethod2.h"
#include "GLMSampler.h"

#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>
#include <sampler/SingletonGraphView.h>
#include <sampler/MutableSampler.h> 

#include <algorithm>
#include <utility>

using std::vector;
using std::string;
using std::list;
using std::set;

namespace jags {
    namespace glm {

	REFactory2::REFactory2(std::string const &name)
	    : _name(name)
	{}
	
	bool REFactory2::checkTau(SingletonGraphView const *tau,
				  GraphView const *glmview) const
	{
	    if (!tau->deterministicChildren().empty()) {
		return false;
	    }
	    
	    vector<StochasticNode *> const &eps = tau->stochasticChildren();
	    for (unsigned int i = 0; i < eps.size(); ++i) {
		if (isObserved(eps[i])) {
		    return false;
		}
		if (isBounded(eps[i])) {
		    return false; 
		}
		if (eps[i]->distribution()->name() != "dnorm" &&
		    eps[i]->distribution()->name() != "dmnorm") {
		    return false;
		}

		Node const *mu_tau = eps[i]->parents()[1];
		if (mu_tau != tau->node()) {
		    return false;
		}
		if (tau->isDependent(eps[i]->parents()[0])) {
		    return false; //mean parameter depends on snode
		}
	    }

	    //Check that all stochastic children of tau are included in
	    //the linear predictor of the glm
	    if (glmview->nodes().size() < eps.size()) {
		return false;
	    }
	    
	    set<StochasticNode*> lpset;
	    lpset.insert(glmview->nodes().begin(), glmview->nodes().end());
	    for (unsigned int i = 0; i < eps.size(); ++i) {
		if (lpset.count(eps[i]) == 0) {
		    return false;
		}
	    }
	    
	    return true; //We made it!
	}
	
	REFactory2::~REFactory2()
	{}
    
	Sampler * REFactory2::makeSampler(
	    list<StochasticNode*> const &free_nodes,
	    set<StochasticNode*> &used_nodes,
	    GLMSampler const *glmsampler, Graph const &graph) const
	{
	    SingletonGraphView *tau = 0;
	    for (list<StochasticNode*>::const_iterator p = free_nodes.begin();
		 p != free_nodes.end(); ++p)
	    {
		if (used_nodes.count(*p)) continue;
		
		if (canSample(*p)) {
		    tau = new SingletonGraphView(*p, graph);
		    if (checkTau(tau, glmsampler->_view)) {
			break;
		    }
		    else {
			delete tau; tau = 0;
		    }
		}
	    }
	    
	    /* Create a single GraphView containing all sampled nodes
	       (from tau and from eps). This is required by the Sampler
	       class. Note that this is a multilevel GraphView.

	    vector<StochasticNode*> snodes = eps->nodes();
	    snodes.push_back(tau->node());
	    GraphView *view = new GraphView(snodes, graph, true);
	    */

	    if (tau) {
		unsigned int nchain = glmsampler->_methods.size();
		vector<MutableSampleMethod*> methods(nchain);
		for (unsigned int i = 0; i < nchain; ++i) {
		    methods[i] = newMethod(tau, glmsampler->_methods[i]);
		}
		used_nodes.insert(tau->node());
		return new MutableSampler(tau, methods, _name);
	    }
	    return 0;
	}
	
    } // namespace glm
} //namespace jags

