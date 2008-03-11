#include <config.h>

#include "ConjugateDirichlet.h"

#include <rng/RNG.h>
#include <graph/MixtureNode.h>
#include <graph/StochasticNode.h>
#include <graph/NodeError.h>
#include "ConjugateSampler.h"

#include <set>
#include <stdexcept>
#include <vector>
#include <cmath>

#include <JRmath.h>

using std::vector;
using std::set;
using std::sqrt;
using std::logic_error;
using std::string;

bool ConjugateDirichlet::canSample(StochasticNode *snode, Graph const &graph)
{
    if(getDist(snode) != DIRCH)
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
	vector<Node const *> const &param = stoch_nodes[i]->parents();
	if (isBounded(stoch_nodes[i])) {
	    return false; //Truncated
	}
	switch(getDist(stoch_nodes[i])) {
	case CAT:
	    break;
	case MULTI:
	    if (param[1] == snode)
		return false;
	    break;
	default:
	    return false;
	}
    }

    // Check deterministic descendants
    for (unsigned int j = 0; j < dtrm_nodes.size(); ++j) {
	MixtureNode const *mnode = asMixture(dtrm_nodes[j]);
	if (mnode) {
	    // Check that indices do not depend on snode
	    unsigned int nindex = mnode->index_size();
	    for (unsigned int i = 0; i < nindex; ++i) {
		if (paramset.count(mnode->parents()[i]))
		    return false;
	    }
	}
	else {
	    return false;
	}
    }
  
    return true;
}

static bool allzero(double const *x, long length)
{
    for (long i = 0; i < length; ++i) {
	if (x[i])
	    return false;
    }
    return true;
}

void ConjugateDirichlet::initialize(ConjugateSampler *sampler,
				    Graph const &graph)
{
}

void ConjugateDirichlet::update(ConjugateSampler *sampler, unsigned int chain, 
				RNG *rng) const
{
    StochasticNode *snode = sampler->node();
    unsigned long size = snode->length();
    double *alpha = new double[size];
    double const *prior = snode->parents()[0]->value(chain);
    for (unsigned long i = 0; i < size; ++i) {
	alpha[i] = prior[i];
    }

    /* Set value to zero:
       This is an illegal value. It is used to find relevant
       children in mixture models.
       FIXME: We can ignore this in non-mixture models
    */
    double *xnew = new double[size];
    for (unsigned long i = 0; i < size; ++i) {
	xnew[i] = 0;
    }
    sampler->setValue(xnew, size, chain);

    vector<StochasticNode const*> const &stoch_children = 
	sampler->stochasticChildren();
    vector<ConjugateDist> const &child_dist = sampler->childDist();
    unsigned int nchildren = stoch_children.size();
    for (unsigned int i = 0; i < nchildren; ++i) {
	StochasticNode const *schild = stoch_children[i];
	long index = 0;
	double const *N = 0;
	if (allzero(schild->parents()[0]->value(chain), 
		    schild->parents()[0]->length())) {
	    switch(child_dist[i]) {
	    case MULTI:
		N = schild->value(chain);
		for (unsigned long j = 0; j < size; ++j) {
		    alpha[j] += N[j];
		}
		break;
	    case CAT:
		index = static_cast<long>(*schild->value(chain) + 1.0E-6);
		alpha[index - 1] += 1;
		break;
	    default:
		throw logic_error("Invalid distribution in Conjugate Dirichlet sampler");
	    }
	}
    }

    /* Check structural zeros */
    for (unsigned long i = 0; i < size; ++i) {
	if (prior[i] == 0 && alpha[i] != 0) {
	    throw NodeError(snode, "Invalid likelihood for Dirichlet distribution with structural zeros");
	}
    }
  
    /* 
       Draw Dirichlet sample by drawing independent gamma random
       variates and then normalizing
    */

    double xsum = 0.0;
    for (unsigned long i = 0; i < size; ++i) {
	if (alpha[i] > 0) {
	    xnew[i] = rgamma(alpha[i], 1, rng);
	    xsum += xnew[i];
	}
	else {
	    xnew[i] = 0;
	}
    }
    for (unsigned long i = 0; i < size; ++i) {
	xnew[i] /= xsum;
    }

    sampler->setValue(xnew, size, chain);

    delete [] xnew;
    delete [] alpha;
}

string ConjugateDirichlet::name() const
{
    return "ConjugateDirichlet";
}
