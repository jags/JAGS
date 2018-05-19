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

static void calBeta(double *betas, SingletonGraphView const *gv,
                    unsigned int chain)
{
    StochasticNode *snode = gv->node();
    double const *xold = snode->value(chain);
    unsigned long nrow = snode->length();

    double *xnew = new double[nrow];
    for (unsigned long i = 0; i < nrow; ++i) {
	xnew[i] = xold[i];
    }

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
	gv->setValue(xnew, nrow, chain);
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
    gv->setValue(xnew, nrow, chain);

    delete [] xnew;
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
    : ConjugateMethod(gv), _betas(0), 
      _length_betas(sumChildrenLength(gv) * gv->length())
{
    if(!gv->deterministicChildren().empty() && checkLinear(gv, true))
    {
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
    #pragma omp critical
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
    double *b = new double[nrow];
    double *A = new double[N];
    for (unsigned long i = 0; i < nrow; ++i) {
	b[i] = 0;
	for (unsigned long i2 = 0; i2 < nrow; ++i2) {
	    b[i] += priorprec[i * nrow + i2] * (priormean[i2] - xold[i2]);
	}
    }
    for (unsigned long i = 0; i < N; ++i) {
	A[i] = priorprec[i];
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

	double *delta = new double[nrow]; 

	int Ni = asInteger(N);
	int ni = asInteger(nrow);

	for (unsigned long j = 0; j < nchildren; ++j) {
	    double const *Y = stoch_children[j]->value(chain);
	    double const *tau = stoch_children[j]->parents()[1]->value(chain);
	    double alpha = 1;

	    F77_DAXPY (&Ni, &alpha, tau, &i1, A, &i1);
	    for (unsigned long i = 0; i < nrow; ++i) {
		delta[i] = Y[i] - xold[i];
	    }
	    F77_DGEMV ("N", &ni, &ni, &alpha, tau, &ni, delta, &i1,
		       &d1, b, &i1);
	}

	delete [] delta;
	
    }
    else {
	
	bool temp_beta = (_betas == 0);
        double *betas = 0;
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
	double *C = new double[nrow * max_nrow_child];
	double *delta = new double[max_nrow_child];
	
	/* Now add the contribution of each term to A, b 
	   
	   b += N_j * beta_j %*% tau_j (Y_j - mu_j)
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

		double alpha = tau[0];
		F77_DSYR("L", &ni, &alpha, beta_j, &i1, A, &ni);
		alpha *= (Y[0] - mu[0]);
		F77_DAXPY(&ni, &alpha, beta_j, &i1, b, &i1);

	    }
	    else {

		double alpha = 1;
		int nc = asInteger(nrow_child);
			    

		F77_DSYMM("R", "L", &ni, &nc, &alpha, tau,
                          &nc, beta_j, &ni, &zero, C, &ni);

		for (unsigned int i = 0; i < nrow_child; ++i) {
		    delta[i] = Y[i] - mu[i];
		}
		F77_DGEMV("N", &ni, &nc, &d1, C, &ni,
			  delta, &i1, &d1, b, &i1);
		F77_DGEMM("N","T", &ni, &ni, &nc,
			  &d1, C, &ni, beta_j, &ni, &d1, A, &ni);
	    }
	       
	    beta_j += nrow_child * nrow;
	}

	delete [] C;
	delete [] delta;

	if (temp_beta) {
	    delete [] betas;
	}
    }


    /* 
       Solve the equation A %*% x = b to get the posterior mean.
       We have to take a copy of A as it is overwritten during
       the call to DPOSV. The result is stored in b
    */
    double * Acopy = new double[N];
    for (int i = 0; i < N; ++i) {
	Acopy[i] = A[i];
    }
    int one = 1;
    int info;
    int ni = asInteger(nrow);
    
    F77_DPOSV ("L", &ni, &one, Acopy, &ni, b, &ni, &info);
    if (info != 0) {
	delete [] Acopy;
	delete [] A;
	delete [] b;
	throwNodeError(snode,
		       "unable to solve linear equations in ConjugateMNormal");
    }

    //Shift origin back to original scale
    for (unsigned long i = 0; i < nrow; ++i) {
	b[i] += xold[i];
    }
    double *xnew = new double[nrow];
    //NB. This uses the lower triangle of A
    DMNorm::randomsample(xnew, b, A, true, nrow, rng);
    _gv->setValue(xnew, nrow, chain);

    delete [] A;
    delete [] Acopy;
    delete [] b;
    delete [] xnew;
    }
}

}}
