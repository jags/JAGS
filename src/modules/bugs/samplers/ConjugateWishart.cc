#include <config.h>

#include "ConjugateWishart.h"
#include "ConjugateSampler.h"
#include "DWish.h"

#include <rng/RNG.h>
#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>

#include <set>
#include <stdexcept>
#include <vector>
#include <cmath>

#include <JRmath.h>

using std::vector;
using std::set;
using std::sqrt;
using std::invalid_argument;
using std::logic_error;
using std::string;

bool ConjugateWishart::canSample(StochasticNode *snode, Graph const &graph)
{
    if (getDist(snode) != WISH)
	return false;

    if (isBounded(snode))
	return false;
  
    vector<StochasticNode const*> stoch_nodes;
    vector<Node*> dtrm_nodes;
    Sampler::classifyChildren(vector<StochasticNode*>(1,snode), 
			      graph, stoch_nodes, dtrm_nodes);
    /* 
       Create a set of nodes containing snode and its deterministic
       descendants for the checks below.
    */
    set<Node const *> paramset;
    paramset.insert(snode);
    for (unsigned int j = 0; j < dtrm_nodes.size(); ++j) {
	paramset.insert(dtrm_nodes[j]);
    }

    // Check stochastic children
    for (unsigned int i = 0; i < stoch_nodes.size(); ++i) {
	vector<Node const*> const &param = stoch_nodes[i]->parents();
	if (isBounded(stoch_nodes[i])) {
	    return false; //Bounded
	}
	switch(getDist(stoch_nodes[i])) {
	case MNORM:
	    if (paramset.count(param[0])) {
		return false; //mean parameter depends on snode
	    }
	    break;
	default:
	    return false;
	}
    }
  
    // Only direct children are allowed
    if (!dtrm_nodes.empty()) {
	return false;
    }

    return true; //We made it!
}

void ConjugateWishart::initialize(ConjugateSampler *sampler,
				  Graph const &graph)
{
}

void ConjugateWishart::update(ConjugateSampler *sampler, unsigned int chain, 
			      RNG *rng) const
{
    vector<StochasticNode const*> const &stoch_children = 
	sampler->stochasticChildren();
    unsigned int nchildren = stoch_children.size();

    vector<Node const *> const &param = sampler->node()->parents();  

    double k = *param[1]->value(chain);
    double const *Rprior = param[0]->value(chain);
    int nrow = param[0]->dim()[0];

    int N = nrow * nrow;
    double *R = new double[N]; 
    for (int i = 0; i < N; ++i) {
	R[i] = Rprior[i];
    }

    vector<ConjugateDist> const &child_dist = sampler->childDist();
    double *delta = new double[nrow];
    for (unsigned int i = 0; i < nchildren; ++i) {
	StochasticNode const *schild = stoch_children[i];
	if (child_dist[i] != MNORM)
	    throw logic_error("Invalid distribution in Conjugate Wishart sampler");
	vector<Node const *> const &cparam = schild->parents();
    
	double const *Y = schild->value(chain);
	double const *mu = cparam[0]->value(chain);

	for (int j = 0; j < nrow; j++) {
	    delta[j] = Y[j] - mu[j];
	}
	for (int j = 0; j < nrow; j++) {
	    for (int l = 0; l < nrow; l++) {
		R[j*nrow + l] += delta[j] * delta[l];
	    }
	}
	k += 1;
    }
    delete [] delta;

    double *xnew = new double[N];
    DWish::randomSample(xnew, N, R, k, nrow, rng);

    delete [] R;
    sampler->setValue(xnew, N, chain);
    delete [] xnew;
}

string ConjugateWishart::name() const
{
    return "ConjugateWishart";
}
