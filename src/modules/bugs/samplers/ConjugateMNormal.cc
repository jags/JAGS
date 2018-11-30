#include <config.h>
#include <rng/RNG.h>
#include <graph/AggNode.h>
#include <graph/MixtureNode.h>
#include <graph/NodeError.h>
#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <sampler/Linear.h>
#include <sampler/SingletonGraphView.h>
#include <module/ModuleError.h>
#include <util/integer.h>

#include "lapack.h"

#include <set>
#include <vector>
#include <cmath>
#include <string>

#include "ConjugateMNormal.h"
#include "DMNorm.h"

#include <JRmath.h>

using std::string;
using std::vector;
using std::set;
using std::sqrt;
using std::string;

namespace jags {
namespace bugs {

    /* Sample partially observed multivariate normal */
    int randomsample(double *x, double const *b, double const *A, unsigned long nrow,
		     vector<bool> const &observed, unsigned long nobs, RNG *rng)
    {
	unsigned long nfree = nrow - nobs;
	unsigned long N = nfree*nfree;

	// Copy Af = A[f,f], bf = b[f] where f represents indices of
	// free (unobserved) elements
	vector<double> Af(N), bf(nfree);
	if (nobs == 0) {
	    //Completely free
	    copy(b, b + nrow, bf.begin());
	    copy(A, A + N, Af.begin());
	}
	else {
	    //Partly observed, use boolean vector "observed" to find
	    //free elements
	    unsigned long p = 0;
	    for (unsigned long i = 0; i < nrow; ++i) {
		if (!observed[i]) {
		    unsigned long q = 0;
		    for (unsigned long j = 0; j < nrow; ++j) {
			if (!observed[j]) {
			    Af[p * nfree + q++] = A[i * nrow + j];
			}
		    }
		    bf[p++] = b[i];
		}
	    }
	}

	int one = 1;
	int info = 0;
	int nf = asInteger(nfree);

	//solve Af %*% x = bf to get posterior mean. The solution will
	//be in bf after the call to dpotrs.
//	F77_DPOTRS("L", &nf, &one, Af.data(), &nf, bf.data(), &nf, &info);
	if (info != 0) return info;
	
	//Af now holds the Cholesky factorization of Af. Use
	//this to generate a multivariate normal random vector with
	//mean 0 and precision Af.
	vector<double> eps(nfree);
	for (unsigned int i = 0; i < nfree; ++i) {
	    eps[i] = rng->normal();
	}
//	F77_DTRSV("L", "T", "N", &nfree, Af.data(), eps.data(), &one);
		  
	// Copy back sampled values
	unsigned long p = 0;
	for (unsigned int i = 0; i < nrow; ++i) {
	    if (!observed[i]) {
		x[i] += bf[p] + eps[p];
		++p;
	    }
	}
	return info;
    }
    
static void calBeta(double *betas, SingletonGraphView const *gv,
                    unsigned int chain)
{
    StochasticNode *snode = gv->node();
    double const *xold = snode->value(chain);
    unsigned long nrow = snode->length();

    vector<double> xnew(nrow);
    copy(xold, xold + nrow, xnew.begin());

    vector<StochasticNode *> const &stoch_children = 
        gv->stochasticChildren();

    unsigned long nchildren = stoch_children.size();
    double *beta_j = betas;
    for (unsigned long j = 0; j < nchildren; ++j) {
	StochasticNode const *schild = stoch_children[j];
	double const *mu = schild->parents()[0]->value(chain);
	unsigned long nrow_child = schild->length();
	for (unsigned long k = 0; k < nrow_child; ++k) {
	    for (unsigned long i = 0; i < nrow; ++i) {
		beta_j[nrow * k + i] = -mu[k];
	    }
	}
	beta_j += nrow_child * nrow;
    }

    for (unsigned int i = 0; i < nrow; ++i) {
	xnew[i] += 1;
	gv->setValue(xnew, chain);
	beta_j = betas;
	for (unsigned long j = 0; j < nchildren; ++j) {
	    StochasticNode const *schild = stoch_children[j];
	    double const *mu = schild->parents()[0]->value(chain);
	    unsigned long nrow_child = schild->length();
	    for (unsigned long k = 0; k < nrow_child; ++k) {
		beta_j[nrow * k + i] += mu[k];
	    }
	    beta_j += nrow_child * nrow;
	}
	xnew[i] -= 1;
    }
    gv->setValue(xnew, chain);
}

static unsigned int sumChildrenLength(SingletonGraphView const *gv)
{
    vector<StochasticNode *> const &children = 
	gv->stochasticChildren(); 

    unsigned int N = 0;
    for (unsigned int i = 0; i < children.size(); ++i) {
	N += children[i]->length();
    }
    return N;
}

ConjugateMNormal::ConjugateMNormal(SingletonGraphView const *gv)
    : ConjugateMethod(gv), _betas(nullptr), 
      _length_betas(sumChildrenLength(gv) * gv->length())
{
    if(checkLinear(gv, true)) {
	_betas = new double[_length_betas];
	calBeta(_betas, gv, 0);
    }
}

ConjugateMNormal::~ConjugateMNormal()
{
    delete [] _betas;
}

bool ConjugateMNormal::canSample(StochasticNode *snode, Graph const &graph)
{
    if (getDist(snode) != MNORM)
	return false;
  
    if (isBounded(snode))
	return false;

    SingletonGraphView gv(snode, graph);
    vector<StochasticNode *> const &schild = gv.stochasticChildren();

    // Check stochastic children
    for (unsigned int i = 0; i < schild.size(); ++i) {
	if (getDist(schild[i]) != MNORM && getDist(schild[i]) != NORM) {
	    return false; //Not normal or multivariate normal
	}
	if (isBounded(schild[i])) {
	    return false;
	}
	if (gv.isDependent(schild[i]->parents()[1])) {
	    return false; //Precision depends on snode
	}
    }

    // Check linearity of deterministic descendants
    if (!checkLinear(&gv, false))
	return false;

    return true; //We made it!
}

void ConjugateMNormal::update(unsigned int chain, RNG *rng) const
{
    vector<StochasticNode *> const &stoch_children = 
          _gv->stochasticChildren();
    unsigned long nchildren = stoch_children.size();
    
    StochasticNode *snode = _gv->node();
    double const *xold = snode->value(chain);
    double const *priormean = snode->parents()[0]->value(chain); 
    double const *priorprec = snode->parents()[1]->value(chain);
    unsigned long nrow = snode->length();
    /* 
       The log of the full conditional density takes the form
       -1/2(t(x) %*% A %*% x - 2 * b %*% x)

       For computational convenience, we reset the origin to xold,
       the current value of the node.
    */
    unsigned long N = nrow * nrow;
    vector<double> A(N);
    copy(priorprec, priorprec + N, A.begin());

    vector<double> b(nrow, 0);
    for (unsigned long i = 0; i < nrow; ++i) {
	for (unsigned long j = 0; j < nrow; ++j) {
	    b[i] += priorprec[i * nrow + j] * (priormean[j] - xold[j]);
	}
    }
    
    /* FORTRAN routines are all call-by-reference, so we need to create
     * these constants */
    double zero = 0;
    double d1 = 1;
    int i1 = 1;
    
    if (_gv->deterministicChildren().empty()) {
      
	// This can only happen if the stochastic children are all
	// multivariate normal with the same number of rows and 
	// columns. We know alpha = 0, beta = I.

	vector<double> delta(nrow);

	int Ni = asInteger(N);
	int ni = asInteger(nrow);

	for (unsigned long j = 0; j < nchildren; ++j) {
	    double const *Y = stoch_children[j]->value(chain);
	    double const *tau = stoch_children[j]->parents()[1]->value(chain);
	    double alpha = 1;

	    F77_DAXPY (&Ni, &alpha, tau, &i1, A.data(), &i1);
	    for (unsigned long i = 0; i < nrow; ++i) {
		delta[i] = Y[i] - xold[i];
	    }
	    F77_DGEMV ("N", &ni, &ni, &alpha, tau, &ni, delta.data(), &i1,
		       &d1, b.data(), &i1);
	}
    }
    else {
	
	bool temp_beta = (_betas == nullptr);
        double *betas = nullptr;
	if (temp_beta) {
	    betas = new double[_length_betas];
	    calBeta(betas, _gv, chain);
	}
        else {
            betas = _betas;
        }

	//Calculate largest possible size of working matrix C
	unsigned long max_nrow_child = 0;
	for (unsigned long j = 0; j < nchildren; ++j) {
	    unsigned long nrow_j = stoch_children[j]->length();
	    if (nrow_j > max_nrow_child) max_nrow_child = nrow_j;
	}
	vector<double> C(nrow * max_nrow_child);
	vector<double> delta(max_nrow_child);
	
	/* Now add the contribution of each term to A, b 
	   
	   b += N_j * beta_j %*% tau_j %*% (Y_j - mu_j)
	   A += N_j * beta_j %*% tau_j %*% t(beta_j)

	   where 
	   - N_j is the frequency weight of child j
	   - beta_j is a matrix of linear coefficients
	   - tau_j is the variance-covariance matrix of child j
	   - mu_j is the mean of child j
	   - Y_j is the value of child j
	   
	   We make use of BLAS routines for efficiency.

	 */
	double const *beta_j = betas;
	for (unsigned long j = 0; j < nchildren; ++j) {
	    
	    StochasticNode const *schild = stoch_children[j];
	    double const *Y = schild->value(chain);
	    double const *mu = schild->parents()[0]->value(chain);
	    double const *tau = schild->parents()[1]->value(chain);
	    unsigned long nrow_child = schild->length();

	    int ni = asInteger(nrow);
	    
	    if (nrow_child == 1) {
		// Scalar children: normal
		double alpha = tau[0];
		F77_DSYR("L", &ni, &alpha, beta_j, &i1, A.data(), &ni);
		alpha *= (Y[0] - mu[0]);
		F77_DAXPY(&ni, &alpha, beta_j, &i1, b.data(), &i1);
	    }
	    else {
		// Vector children: multivariate normal
		double alpha = 1;
		int nc = asInteger(nrow_child);

		F77_DSYMM("R", "L", &ni, &nc, &alpha, tau,
                          &nc, beta_j, &ni, &zero, C.data(), &ni);

		for (unsigned int i = 0; i < nrow_child; ++i) {
		    delta[i] = Y[i] - mu[i];
		}
		F77_DGEMV("N", &ni, &nc, &d1, C.data(), &ni,
			  delta.data(), &i1, &d1, b.data(), &i1);
		F77_DGEMM("N","T", &ni, &ni, &nc,
			  &d1, C.data(), &ni, beta_j, &ni, &d1, A.data(), &ni);
	    }
	       
	    beta_j += nrow_child * nrow;
	}

	if (temp_beta) {
	    delete [] betas;
	}
    }


    /* 
       Solve the equation A %*% x = b to get the posterior mean.
       We have to take a copy of A as it is overwritten during
       the call to DPOSV. The result is stored in b
    */
    vector<double> F(A);

    int one = 1;
    int info;
    int ni = asInteger(nrow);
    
    F77_DPOSV ("L", &ni, &one, F.data(), &ni, b.data(), &ni, &info);
    if (info != 0) {
	throwNodeError(snode,
		       "unable to solve linear equations in ConjugateMNormal");
    }

    //Shift origin back to original scale
    for (unsigned long i = 0; i < nrow; ++i) {
	b[i] += xold[i];
    }
    vector<double> xnew(nrow);
    //NB. This uses the lower triangle of A
    DMNorm::randomsample(xnew.data(), b.data(), A.data(), true, nrow, rng);
    _gv->setValue(xnew, chain);
}

}}
