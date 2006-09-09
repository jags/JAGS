#include <config.h>
#include <rng/RNG.h>
#include <graph/AggNode.h>
#include <graph/MixtureNode.h>
#include <graph/NodeError.h>
#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include "lapack.h"

#include <set>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <string>

#include "ConjugateMNormal.h"
#include "ConjugateSampler.h"
#include "DMNorm.h"

#include <Rmath.h>

using std::string;
using std::vector;
using std::set;
using std::sqrt;
using std::invalid_argument;

ConjugateMNormal::ConjugateMNormal(StochasticNode *snode, 
				   Graph const &graph, unsigned int chain)
    : ConjugateSampler(snode, graph, chain), _betas(0), _length_betas(0)
{
    if (!canSample(snode, graph)) {
	string msg("Can't construct ConjugateMNormal sampler");
	throw invalid_argument(msg);
    }

    if (!deterministicChildren().empty()) {

	vector<StochasticNode const *> const &children = stochasticChildren();

	unsigned int N = 0;
	for (unsigned int i = 0; i < children.size(); ++i) {
	    N += children[i]->length();
	}
	_length_betas = N * snode->length();

	// Check for constant linear terms
	bool fixed_beta = true;
	set<Node*> paramset;
	vector<DeterministicNode*> const &dtrm = deterministicChildren();
	paramset.insert(snode);
	for (unsigned int j = 0; j < dtrm.size(); ++j) {
	    paramset.insert(dtrm[j]);
	}
	for (unsigned int j = 0; j < dtrm.size(); ++j) {
	    if (!dtrm[j]->isLinear(paramset, true)) {
		fixed_beta = false;
		break;
	    }
	}
	
	if (fixed_beta) {
	    _betas = new double[_length_betas];
	    calBeta();
	}
    }
}

ConjugateMNormal::~ConjugateMNormal()
{
    delete [] _betas;
}


void ConjugateMNormal::calBeta()
{
    double const *xold = node()->value(chain());
    unsigned int nrow = node()->length();

    double *xnew = new double[nrow];
    for (int i = 0; i < nrow; ++i) {
	xnew[i] = xold[i];
    }

    vector<StochasticNode const*> const &stoch_children = stochasticChildren();

    unsigned long nchildren = stoch_children.size();
    double *beta_j = _betas;
    for (unsigned int j = 0; j < nchildren; ++j) {
	StochasticNode const *snode = stoch_children[j];
	double const *mu = snode->parents()[0]->value(chain());
	unsigned int nrow_child = snode->length();
	for (int k = 0; k < nrow_child; ++k) {
	    for (unsigned int i = 0; i < nrow; ++i) {
		beta_j[nrow * k + i] = -mu[k];
	    }
	}
	beta_j += nrow_child * nrow;
    }

    for (int i = 0; i < nrow; ++i) {
	xnew[i] += 1;
	setValue(xnew, nrow);
	beta_j = _betas;
	for (unsigned int j = 0; j < nchildren; ++j) {
	    StochasticNode const *snode = stoch_children[j];
	    double const *mu = snode->parents()[0]->value(chain());
	    unsigned int nrow_child = snode->length();
	    for (int k = 0; k < nrow_child; ++k) {
		beta_j[nrow * k + i] += mu[k];
	    }
	    beta_j += nrow_child * nrow;
	}
	xnew[i] -= 1;
    }
    setValue(xnew, nrow);

    delete [] xnew;
}


bool ConjugateMNormal::canSample(StochasticNode *snode, Graph const &graph)
{
    if (getDist(snode) != MNORM)
	return false;
  
    if (snode->isBounded())
	return false;

    vector<StochasticNode const*> stoch_nodes;
    vector<DeterministicNode*> dtrm_nodes, extra_nodes;
    classifyChildren(vector<StochasticNode*>(1,snode), 
		     graph, stoch_nodes, dtrm_nodes);
    /* 
       Create a set of nodes containing snode and its deterministic
       descendants for the checks below.
    */
    set<Node *> paramset;
    paramset.insert(snode);
    for (unsigned int j = 0; j < dtrm_nodes.size(); ++j) {
	paramset.insert(dtrm_nodes[j]);
    }

    // Check stochastic children
    for (unsigned int i = 0; i < stoch_nodes.size(); ++i) {
	if (getDist(stoch_nodes[i]) != MNORM &&
	    getDist(stoch_nodes[i]) != NORM) {
	    return false; //Not normal or multivariate normal
	}
	if (stoch_nodes[i]->isBounded()) {
	    return false;
	}
	vector<Node *> const &param = stoch_nodes[i]->parents();
	if (paramset.count(param[1])) {
	    return false; //Precision depends on snode
	}
    }

    // Check deterministic descendants
    for (unsigned int j = 0; j < dtrm_nodes.size(); ++j) {
	if (!dtrm_nodes[j]->isLinear(paramset, false)) {
	    return false;
	}
    }
  
    return true; //We made it!
}

void ConjugateMNormal::update(RNG *rng)
{
    vector<StochasticNode const*> const &stoch_children = stochasticChildren();
    unsigned int nchildren = stoch_children.size();
    
    double const *xold = node()->value(chain());
    double const *priormean = node()->parents()[0]->value(chain()); 
    double const *priorprec = node()->parents()[1]->value(chain());
    int nrow = node()->length();
    /* 
       The log of the full conditional density takes the form
       -1/2(t(x) %*% A %*% x - 2 * b %*% x)

       For computational convenience, we reset the origin to xold,
       the current value of the node.
    */
    int N = nrow * nrow;
    double *b = new double[nrow];
    double *A = new double[N];
    for (int i = 0; i < nrow; ++i) {
	b[i] = 0;
	for (int i2 = 0; i2 < nrow; ++i2) {
	    b[i] += priorprec[i * nrow + i2] * (priormean[i2] - xold[i2]);
	}
    }
    for (int i = 0; i < N; ++i) {
	A[i] = priorprec[i];
    }
    
    if (deterministicChildren().empty()) {
      
	// This can only happen if the stochastic children are all
	// multivariate normal with the same number of rows and 
	// columns normal. We know alpha = 0, beta = I.
	
	for (unsigned int j = 0; j < nchildren; ++j) {
	    double const *Y = stoch_children[j]->value(chain());
	    double const *tau = stoch_children[j]->parents()[1]->value(chain());
	    for (int i = 0; i < nrow; ++i) {
		for (int i2 = 0; i2 < nrow; ++i2) {
		    A[i * nrow + i2] += tau[i * nrow + i2];
		    b[i] += (Y[i2] - xold[i2]) * tau[nrow * i + i2];
		}
	    }
	}
	
    }
    else {
	
	bool temp_beta = (_betas == 0);
	if (temp_beta) {
	    _betas = new double[_length_betas];
	    calBeta();
	}
	
	/* Now add the contribution of each term to A, b */
	double const *beta_j = _betas;
	for (unsigned int j = 0; j < nchildren; ++j) {
	    
	    Node const *snode = stoch_children[j];
	    double const *Y = snode->value(chain());
	    double const *mu = snode->parents()[0]->value(chain());
	    double const *tau = snode->parents()[1]->value(chain());
	    int nrow_child = snode->length();
	  
	    for (int i = 0; i < nrow; ++i) {
		for (int k = 0; k < nrow_child; ++k) {
		  
		    double beta_tau = 0;
		    for (int k2 = 0; k2 < nrow_child; ++k2) {
			beta_tau += 
			    beta_j[nrow * k2 + i] * tau[nrow_child * k2 + k];
		    }
		  
		    for (int i2 = 0; i2 < nrow; ++i2) {
			A[i * nrow + i2] += beta_tau * beta_j[nrow * k + i2];
		    }
		  
		    b[i] += (Y[k] - mu[k]) * beta_tau;
		  
		}
	    }
	  
	    beta_j += nrow_child * nrow;
	}

	if (temp_beta) {
	    delete _betas;
	    _betas = 0;
	}
    }


    /* 
       Solve the equation A %*% x = b to get the posterior mean.
       We have to take a copy of A as it is overwritten during
       the call to DSYSV. The result is stored in b
    */
    double * Acopy = new double[N];
    for (int i = 0; i < N; ++i) {
	Acopy[i] = A[i];
    }
    int *ipiv = new int[nrow];
    int one = 1;
    double worktest;
    int lwork = -1;
    int info;
    F77_DSYSV ("L", &nrow, &one, Acopy, &nrow, ipiv, b, &nrow, &worktest, 
	       &lwork, &info);
    if (info != 0) {
	throw NodeError(node(),
			"unable to solve linear equations in Conjugate mnorm sampler");
  }
    lwork = static_cast<int>(worktest) + 1;
    double * work = new double[lwork];
    F77_DSYSV ("L", &nrow, &one, Acopy, &nrow, ipiv, b, &nrow, work, &lwork,
	       &info);
    if (info != 0) {
	throw NodeError(node(),
			"unable to solve linear equations in Conjugate MNorm sampler");
    }
    delete [] work;
    delete [] Acopy;
    delete [] ipiv;


    //Shift origin back to original scale
    for (unsigned int i = 0; i < nrow; ++i) {
	b[i] += xold[i];
    }
    double *xnew = new double[nrow];
    DMNorm::randomsample(xnew, b, A, nrow, rng);
    setValue(xnew, nrow);
    
    delete [] b;
    delete [] A;
    delete [] xnew;
}
