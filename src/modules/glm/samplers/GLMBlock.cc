#include <config.h>

#include "GLMBlock.h"
#include "Outcome.h"

#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>

#include <cmath>

using std::vector;
using std::sqrt;

extern cholmod_common *glm_wk;

namespace jags {

namespace glm {

    GLMBlock::GLMBlock(GraphView const *view, 
			 vector<SingletonGraphView const *> const &sub_views,
			 vector<Outcome *> const &outcomes,
			 unsigned int chain)
	: GLMMethod(view, sub_views, outcomes, chain)
    {
	calDesign();
	symbolic();
    }

    void GLMBlock::update(RNG *rng) 
    {
	//   The log of the full conditional density takes the form
	//   -(t(x) %*% A %*% x - 2 * b %*% x)/2
	//   where A is the posterior precision and the mean mu solves
	//   A %*% mu = b

	//   For computational convenience we take xold, the current value
	//   of the sampled nodes, as the origin

	// Update outcomes
	for (vector<Outcome*>::const_iterator p = _outcomes.begin();
	     p != _outcomes.end(); ++p)
	{
	    (*p)->update(rng);
	}
	
	double *b = 0;
	cholmod_sparse *A = 0;
	calCoef(b, A);
	
	// Get LDL' decomposition of posterior precision
	A->stype = -1;
        #pragma omp critical
	{
	    int ok = cholmod_factorize(A, _factor, glm_wk);
	    cholmod_free_sparse(&A, glm_wk);
	    if (!ok) {
		throwRuntimeError("Cholesky decomposition failure in GLMBlock");
	    }
	}

	// Use the LDL' decomposition to generate a new sample
	// with mean mu such that A %*% mu = b and precision A. 
	
	unsigned int nrow = _view->length();
	cholmod_dense *w = 0;
	#pragma omp critical
	w = cholmod_allocate_dense(nrow, 1, nrow, CHOLMOD_REAL, glm_wk);

	// Permute RHS
	double *wx = static_cast<double*>(w->x);
	int *perm = static_cast<int*>(_factor->Perm);
	for (unsigned int i = 0; i < nrow; ++i) {
	    wx[i] = b[perm[i]];
	}

	cholmod_dense *u1 = 0;
        #pragma omp critical
	u1 = cholmod_solve(CHOLMOD_L, _factor, w, glm_wk);
	updateAuxiliary(u1, _factor, rng);

	double *u1x = static_cast<double*>(u1->x);
	if (_factor->is_ll) {
	    // LL' decomposition
	    for (unsigned int r = 0; r < nrow; ++r) {
		u1x[r] += rng->normal();
	    }
	}
	else {
	    // LDL' decomposition. The diagonal D matrix is stored
	    // as the diagonal of _factor
	    int *fp = static_cast<int*>(_factor->p);
		double *fx = static_cast<double*>(_factor->x);
		for (unsigned int r = 0; r < nrow; ++r) {
		    u1x[r] += rng->normal() * sqrt(fx[fp[r]]);
		}
	}

	cholmod_dense *u2 = 0;
        #pragma omp critical
	u2 = cholmod_solve(CHOLMOD_DLt, _factor, u1, glm_wk);

	// Permute solution
	double *u2x = static_cast<double*>(u2->x);
	for (unsigned int i = 0; i < nrow; ++i) {
	    b[perm[i]] = u2x[i];
	}

        #pragma omp critical
	{
            cholmod_free_dense(&w, glm_wk);
            cholmod_free_dense(&u1, glm_wk);
            cholmod_free_dense(&u2, glm_wk);
        }
	//Shift origin back to original scale
	int r = 0;
	for (vector<StochasticNode*>::const_iterator p = 
		 _view->nodes().begin();  p != _view->nodes().end(); ++p)
	{
	    unsigned int length = (*p)->length();
	    double const *xold = (*p)->value(_chain);
	    for (unsigned int i = 0; i < length; ++i, ++r) {
		b[r] += xold[i];
	    }
	}

	_view->setValue(b, nrow, _chain);
	delete [] b;
    }

    void GLMBlock::updateAuxiliary(cholmod_dense *b, cholmod_factor *N,
				   RNG *rng)
    {
    }

}}
