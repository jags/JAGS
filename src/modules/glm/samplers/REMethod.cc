#include <config.h>

#include "REMethod.h"
#include "Outcome.h"

#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <JRmath.h>

#include <cmath>

using std::vector;
using std::sqrt;

extern cholmod_common *glm_wk;

namespace jags {
    namespace glm {

	REMethod::REMethod(SingletonGraphView const *tau,
			   GraphView const *eps,
			   vector<SingletonGraphView const *> const &sub_eps,
			   vector<Outcome *> const &outcomes,
			   unsigned int chain)
	    : GLMMethod(eps, sub_eps, outcomes, chain), _tau(tau)
	{
	    calDesign();
	    symbolic();
	}

	//FIXME: This is largely copy-pasted from GLMBlock. Surely no need
	//to reproduce it here?
	void REMethod::updateEps(RNG *rng) 
	{
	    //   The log of the full conditional density takes the form
	    //   -(t(x) %*% A %*% x - 2 * b %*% x)/2
	    //   where A is the posterior precision and the mean mu solves
	    //   A %*% mu = b

	    //   For computational convenience we take xold, the
	    //   current value of the sampled nodes, as the origin

	
	    double *b = 0;
	    cholmod_sparse *A = 0;
	    calCoef(b, A);
	
	    // Get LDL' decomposition of posterior precision
	    A->stype = -1;
	    int ok = cholmod_factorize(A, _factor, glm_wk);
	    cholmod_free_sparse(&A, glm_wk);
	    if (!ok) {
		throwRuntimeError("Cholesky decomposition failure in REMethod");
	    }

	    // Use the LDL' decomposition to generate a new sample
	    // with mean mu such that A %*% mu = b and precision A. 
	
	    unsigned int nrow = _view->length();
	    cholmod_dense *w =
		cholmod_allocate_dense(nrow, 1, nrow, CHOLMOD_REAL, glm_wk);

	    // Permute RHS
	    double *wx = static_cast<double*>(w->x);
	    int *perm = static_cast<int*>(_factor->Perm);
	    for (unsigned int i = 0; i < nrow; ++i) {
		wx[i] = b[perm[i]];
	    }

	    cholmod_dense *u1 = cholmod_solve(CHOLMOD_L, _factor, w, glm_wk);
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

	    cholmod_dense *u2 = cholmod_solve(CHOLMOD_DLt, _factor, u1, glm_wk);

	    // Permute solution
	    double *u2x = static_cast<double*>(u2->x);
	    for (unsigned int i = 0; i < nrow; ++i) {
		b[perm[i]] = u2x[i];
	    }

	    cholmod_free_dense(&w, glm_wk);
	    cholmod_free_dense(&u1, glm_wk);
	    cholmod_free_dense(&u2, glm_wk);

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

	void REMethod::update(RNG *rng) {

	    // Update outcomes
	    for (vector<Outcome*>::const_iterator p = _outcomes.begin();
		 p != _outcomes.end(); ++p)
	    {
		(*p)->update(rng);
	    }

	    updateEps(rng);
	    updateTau(rng);
	
	}

    }
}
