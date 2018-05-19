#include <config.h>

#include "HolmesHeld.h"
#include "Outcome.h"

#include <graph/StochasticNode.h>
#include <graph/LinkNode.h>
#include <sampler/GraphView.h>
#include <rng/TruncatedNormal.h>
#include <rng/RNG.h>
#include <module/ModuleError.h>

#include <cmath>

// Regularization penalty for precision 
#define REG_PENALTY 0.001

using std::vector;
using std::string;
using std::sqrt;

extern cholmod_common *glm_wk;

static cholmod_sparse shallow_copy(cholmod_sparse *x, unsigned int c)
{
    //Take a copy of column c of sparse matrix x without allocating
    //any memory. This is computationally cheaper than calling
    //cholmod_submatrix, but potentially dangerous if the copy is
    //passed to a function that tries to modify it. Note that we can
    //only have one shallow copy in use at a time since they all share
    //the same array of column pointers given by the static array p.

    static int p[2]; 
    
    cholmod_sparse xcopy = *x;

    double *xx = static_cast<double*>(x->x);
    int *xp = static_cast<int*>(x->p);
    int *xi = static_cast<int*>(x->i);

    xcopy.ncol = 1;
    int nz = xp[c+1] - xp[c];
    xcopy.nzmax = nz;
    p[0] = 0;
    p[1] = nz; 
    xcopy.p = &p;
    xcopy.i = xi + xp[c];
    xcopy.x = xx + xp[c];

    return xcopy;
}

namespace jags {
namespace glm {

    
    HolmesHeld::HolmesHeld(GraphView const *view,
			   vector<SingletonGraphView const *> const &sub_views,
			   vector<Outcome *> const &outcomes,
			   unsigned int chain)
	: GLMBlock(view, sub_views, outcomes, chain)
    {
    }

    void HolmesHeld::updateAuxiliary(cholmod_dense *W, 
				     cholmod_factor *N, RNG *rng)
    {
	/* 
	   In the parent GLMMethod class, the posterior precision of
           the regression parameters is represented by the matrix "A";
           the posterior mean "mu" solves A %*% mu = b.
	   
           In this call, "N" holds the factorization of P %*% A %*% t(P), 
           where P is a permutation matrix.  The factorization takes
           the form L %*% D %*% t(L), where D is diagonal and L is
           a lower triangular matrix. The parameter "w" solves
           L %*% w = P %*% b
	   
	   IMPORTANT NOTE: mu, b use a parameterization in which the
	   current value of the regression parameters is taken as the
	   origin.  This requires us to adjust the calculations using
	   the current value of the linear predictor (mu_r).
	*/

	vector<StochasticNode *> const &schildren = 
	    _view->stochasticChildren();

	int nrow = schildren.size();

	cholmod_sparse *pt_x = 0;
	#pragma omp critical
	{
	    //Transpose and permute the design matrix
	    cholmod_sparse *t_x = cholmod_transpose(_x, 1, glm_wk);
	    int *fperm = static_cast<int*>(_factor->Perm);
	    pt_x = cholmod_submatrix(t_x, fperm, t_x->nrow,
				     0, -1, 1, 1, glm_wk);
	    cholmod_free_sparse(&t_x, glm_wk);
	}
	
	int ncol = _x->ncol;
	vector<double> d(ncol, 1);
	if (!_factor->is_ll) {
	    // LDL' decomposition. The diagonal D matrix is stored as
	    // the diagonal of _factor
	    int *fp = static_cast<int*>(_factor->p);
	    double *fx = static_cast<double*>(_factor->x);
	    for (int r = 0; r < ncol; ++r) {
		d[r] = fx[fp[r]];
	    }
	}

	double *Wx = static_cast<double *>(W->x);
	
	cholmod_dense *U = 0, *Y = 0, *E = 0;
	cholmod_sparse *uset = 0;

	cholmod_dense *X = 0;
	#pragma omp critical
	X = cholmod_allocate_dense(ncol, 1, ncol, CHOLMOD_REAL, glm_wk);
	double *Xx = static_cast<double*>(X->x);

	for (int r = 0; r < nrow; ++r) {

	    // We need to do this loop with minimal memory
	    // allocation. The cholmod_solve2 function uses a
	    // dense/sparse pair X/xset for input and another
	    // dense/sparse pair U/uset for output: U, uset, and
	    // workspaces Y and E are allocated on first use and then
	    // reused in subsequent calls (See CHOLMOD
	    // documentation). We also use the utility function
	    // shallow_copy to access each column of pt_x without
	    // further memory allocation.
	    
	    if (_outcomes[r]->fixedb()) continue;

	    cholmod_sparse xset = shallow_copy(pt_x, r);
	    double *xx = static_cast<double*>(xset.x);
	    int *xp = static_cast<int*>(xset.p);
	    int *xi = static_cast<int*>(xset.i);

	    //Copy non-zero elements of column r to dense vector X
	    for (int j = 0; j < xp[1]; ++j) {
		int c = xi[j];
		Xx[c] = xx[j];
	    }

	    #pragma omp critical
	    cholmod_solve2(CHOLMOD_L, _factor, X, &xset, &U, &uset, &Y, &E,
			   glm_wk);

	    double mu_r = _outcomes[r]->mean(); // See IMPORTANT NOTE above
	    double tau_r = _outcomes[r]->precision();
	    
	    //Calculate mean and precision of z[r] conditional
	    //on z[s] for s != r
	    double zr_mean = 0;
	    double gr = 0; 

	    int *up = static_cast<int*>(uset->p);
	    int *ui = static_cast<int*>(uset->i);
	    double *Ux = static_cast<double*>(U->x);

	    for (int j = 0; j < up[1]; ++j) {
		int c = ui[j];
		zr_mean  += Ux[c] * Wx[c] / d[c]; 
		gr  += Ux[c] * Ux[c] / d[c];
	    }

	    double Hr = gr * tau_r; // diagonal of hat matrix
	    if (1 - Hr <= 0) {
		continue;
	    }
	    zr_mean -= Hr * (_outcomes[r]->value() - mu_r);
	    zr_mean /= (1 - Hr);
	    double zr_var = gr / (1 - Hr);
	    
	    double zold = _outcomes[r]->value();
	    _outcomes[r]->update(zr_mean + mu_r, zr_var, rng); 
	    
	    //Add new contribution of row r back to b
	    double zdelta = (_outcomes[r]->value() - zold) * tau_r;
	    for (int j = 0; j < up[1]; ++j) {
		int c = ui[j];
		Wx[c] += Ux[c] * zdelta; 
	    }

	}
	    
	//Free workspace

	#pragma omp critical
	{
	    cholmod_free_sparse(&pt_x, glm_wk);
	    cholmod_free_sparse(&uset, glm_wk);
	    
	    cholmod_free_dense(&U, glm_wk);
	    cholmod_free_dense(&Y, glm_wk);
	    cholmod_free_dense(&E, glm_wk);
	    cholmod_free_dense(&X, glm_wk);
	}
    }
    
}}
    
