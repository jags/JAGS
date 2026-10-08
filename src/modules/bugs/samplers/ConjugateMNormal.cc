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
#include <matrix/blas.h>
#include <matrix/lapack.h>

#include <set>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

#include "ConjugateMNormal.h"
#include "DMNorm.h"
#include "LMDesign.h"

#include <JRmath.h>

using std::string;
using std::vector;
using std::set;
using std::sqrt;
using std::string;
using std::count;

#ifdef DEBUG
#include <module/ModuleError.h>
#endif

namespace jags {
namespace bugs {

    /* 
       Sample partially observed multivariate normal.
    */
    static int MNormSample(vector<double> &x,
			   vector<double> &b, vector<double> &A,
			   vector<bool> const &observed, RNG *rng)
    {
	unsigned long nrow = x.size();
	unsigned long nfree = count(observed.begin(), observed.end(), false);

	if (nfree < nrow) {
	    //Pack the leading parts of matrix A and vector b with
	    //elements corresponding to the free elements of x
	    for (unsigned long i = 0, p = 0; i < nrow; ++i) {
		if (!observed[i]) {
		    for (unsigned long j = 0, q = 0; j < nrow; ++j) {
			if (!observed[j]) {
			    A[p * nrow + q++] = A[i * nrow + j];
			}
		    }
		    b[p++] = b[i];
		}
	    }
	}

	int one = 1;
	int info = 0;
	int nf = asInteger(nfree);
	int nr = asInteger(nrow);

	//Solve A %*% x = b to get the posterior mean. The solution
	//will be in b after the call to dposv.
	jags_dposv("L", &nf, &one, A.data(), &nr, b.data(), &nr, &info);
	if (info != 0) return info;
	
	//After dposv, the leading nfree x nfree lower triangle of A
	//holds the Cholesky factorization. Use it to generate a
	//multivariate normal random vector with mean 0 and precision A
	vector<double> eps(nfree);
	for (unsigned long p = 0; p < nfree; ++p) {
	    eps[p] = rng->normal();
	}
	jags_dtrsv("L", "T", "N", &nf, A.data(), &nr, eps.data(), &one);
		  
	// Copy back sampled values
	for (unsigned long i = 0, p = 0; i < nrow; ++i) {
	    if (!observed[i]) {
		x[i] += b[p] + eps[p];
		++p;
	    }
	}
	return info;
    }

#ifdef DEBUG    
static void calBetaNumeric(double *betas, SingletonGraphView const *gv,
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
#endif
    
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
    : ConjugateMethod(gv), _betas(nullptr), _design(nullptr), 
      _length_betas(sumChildrenLength(gv) * gv->length())
{
    _design = new LMDesign(gv);

    if(checkLinear(gv, true)) {
	//Fixed design matrix
	_betas = new double[_length_betas];
	_design->calCoef(_betas, 0, true);
	delete _design;
	_design = nullptr;
	
#ifdef DEBUG
	vector<double> betas_num(_length_betas);
	calBetaNumeric(betas_num.data(), gv, 0);
	for (unsigned int i = 0; i < _length_betas; ++i) {
	    if (fabs(_betas[i] - betas_num[i]) > 1.0E-6) {
		throwLogicError("ConjugateMNormal::Error calculating design matrix");
	    }
	}
#endif
	
    }
}

ConjugateMNormal::~ConjugateMNormal()
{
    delete [] _betas;
    delete _design;
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

	    jags_daxpy (&Ni, &alpha, tau, &i1, A.data(), &i1);
	    for (unsigned long i = 0; i < nrow; ++i) {
		delta[i] = Y[i] - xold[i];
	    }
	    jags_dgemv ("N", &ni, &ni, &alpha, tau, &ni, delta.data(), &i1,
		       &d1, b.data(), &i1);
	}
    }
    else {
	
	bool temp_beta = (_betas == nullptr);
        double *betas = nullptr;
	if (temp_beta) {
	    /* Design matrix is not fixed and therefore different
	     * between chains and iterations */
	    betas = new double[_length_betas];
	    _design->calCoef(betas, chain, true);
#ifdef DEBUG
	    vector<double> betas_num(_length_betas);
	    calBetaNumeric(betas_num.data(), _gv, chain);
	    for (unsigned int i = 0; i < _length_betas; ++i) {
		if (fabs(betas[i] - betas_num[i]) > 1.0E-6) {
		    throwLogicError("ConjugateMNormal::Error calculating design matrix");
		}
	    }
#endif
	}
	else {
	    /* Fixed design matrix calculated in the constructor */
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
		jags_dsyr("L", &ni, &alpha, beta_j, &i1, A.data(), &ni);
		alpha *= (Y[0] - mu[0]);
		jags_daxpy(&ni, &alpha, beta_j, &i1, b.data(), &i1);
	    }
	    else {
		// Vector children: multivariate normal
		double alpha = 1;
		int nc = asInteger(nrow_child);

		jags_dsymm("R", "L", &ni, &nc, &alpha, tau,
			   &nc, beta_j, &ni, &zero, C.data(), &ni);

		for (unsigned int i = 0; i < nrow_child; ++i) {
		    delta[i] = Y[i] - mu[i];
		}
		jags_dgemv("N", &ni, &nc, &d1, C.data(), &ni,
			   delta.data(), &i1, &d1, b.data(), &i1);
		jags_dgemm("N","T", &ni, &ni, &nc,
			   &d1, C.data(), &ni, beta_j, &ni, &d1, A.data(), &ni);
	    }
	       
	    beta_j += nrow_child * nrow;
	}

	if (temp_beta) {
	    delete [] betas;
	}
    }
   
    vector<double> xnew(nrow);
    copy(xold, xold + nrow, xnew.begin());

    // Note that A, b are overwritten here
    MNormSample(xnew, b, A, *snode->observedMask(), rng);
    _gv->setValue(xnew, chain);
}

}}
