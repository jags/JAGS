#include <config.h>

#include "REMethod.h"
#include "Outcome.h"

#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <JRmath.h>

#include <cmath>
#include <algorithm>

using std::vector;
using std::sqrt;
using std::fill;

extern cholmod_common *glm_wk;

namespace jags {
    namespace glm {

	REMethod::REMethod(SingletonGraphView const *tau,
			   GraphView const *eps,
			   vector<SingletonGraphView const *> const &sub_eps,
			   vector<Outcome *> const &outcomes,
			   unsigned int chain)
	    : GLMMethod(eps, sub_eps, outcomes, chain), _tau(tau), _eps(eps)
	{
	    calDesign();
	    symbolic();

	    unsigned int nrow = _outcomes.size();
	    unsigned int ncol = eps->nodes()[0]->length();
	    _z = cholmod_allocate_dense(nrow, ncol, nrow, CHOLMOD_REAL,
					glm_wk);
	}

	REMethod::~REMethod()
	{
	    cholmod_free_dense(&_z, glm_wk);
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

	void REMethod::calDesignSigma()
	{
	    //Sanity checks
	    unsigned int Neps = _eps->nodes().size();
	    if (_x->nrow != _outcomes.size() || _z->nrow != _x->nrow) {
		throwLogicError("Row mismatch in REMethod");
	    }
	    if (_x->ncol != _z->ncol * Neps || _x->ncol != _eps->length()) {
		throwLogicError("Column mismatch in REMethod");
	    }
	    
	    //Set up access to sparse design matrix for eps
	    int *Xp = static_cast<int*>(_x->p);
	    int *Xi = static_cast<int*>(_x->i);
	    double *Xx = static_cast<double*>(_x->x);
	    
	    //Set up access to dense design matrix for sigma
	    double *Zx = static_cast<double*>(_z->x);

	    //Set all elements of _z to zero.
	    fill(Zx, Zx + _z->nzmax, 0);

	    //If there are m columns of _z then _z[,i] is the sum of
	    //every mth column of _x (starting with _x[,i]),
	    //multiplied by the corresponding random effect

	    vector<StochasticNode*> const &eps=_eps->nodes();
	    for (unsigned int i = 0; i < Neps; ++i) {
		double const *eval = eps[i]->value(_chain);
		double const *emean = eps[i]->parents()[0]->value(_chain);
		for (unsigned int zcol = 0; zcol < _z->ncol; ++zcol) {
		    int xcol = i * _z->ncol + zcol;
		    for (int xi = Xp[xcol]; xi < Xp[xcol+1]; ++xi) {
			int row = Xi[xi];
			int zi = _z->nrow * zcol + row;
			Zx[zi] += Xx[xi] * (eval[zcol] - emean[zcol]);
		    }
		}
	    }
	}

	void REMethod::update(RNG *rng) {
	    
	    // Update outcomes
	    for (vector<Outcome*>::const_iterator p = _outcomes.begin();
		 p != _outcomes.end(); ++p)
	    {
		(*p)->update(rng);
	    }
	    
	    updateEps(rng); //Update random effects
	    updateTau(rng); //Sufficient parameterization
	    updateSigma(rng); //Ancillary parameterization
	    updateTau(rng); //Sufficient parameterization
	}

	void REMethod::calCoefSigma(double *A, double *b, double const *sigma0,
				    unsigned int m) const
	{
	    double const *Zx = static_cast<double const *>(_z->x);
	    unsigned int N = _outcomes.size();

	    for (unsigned int i = 0; i < N; ++i) {
		double Y = _outcomes[i]->value();
		double mu = _outcomes[i]->mean();
		double lambda = _outcomes[i]->precision();
		vector<double> X(m);
		for (unsigned int j = 0; j < m; ++j) {
		    X[j] =  Zx[j*N+i]/sigma0[j];
		}
		for (unsigned int j = 0; j < m; ++j) {
		    for (unsigned int k = 0; k < m; ++k) {
			A[j*m + k] += X[j] * X[k] * lambda;
		    }
		    b[j] += (Y - mu) * X[j] * lambda;
		}
	    }
	}
	
	double REMethod::logLikelihoodSigma(double const *sigma,
					    double const *sigma0,
					    unsigned int m) const
	{
	    vector<double> A(m*m);
	    vector<double> b(m);
	    calCoefSigma(&A[0], &b[0], sigma0, m);

	    vector<double> delta(m);
	    for (unsigned int i = 0; i < m; ++i) {
		delta[i] = sigma[i] - sigma0[i];
	    }
	    
	    double loglik = 0;
	    for (unsigned int i = 0; i < m; ++i) {
		loglik += b[i] * delta[i];
		for (unsigned int j = 0; j < m; ++j) {
		    loglik -= delta[i] * A[i*m + j] * delta[j] /2.0;
		}
	    }
	    return loglik;
	}

	void REMethod::rescaleSigma(double const *sigma,
				    double const *sigma0,
				    unsigned int m)
	{
	    vector<double> sigma_ratio(m);
	    for (unsigned int i = 0; i < m; ++i) {
		sigma_ratio[i] = sigma[i]/sigma0[i];
	    }
	    
	    vector<StochasticNode *> const &eps = _eps->nodes();
	    vector<double> eval(_eps->length());
	    for (unsigned int i = 0; i < eps.size(); ++i) {
		double const *Y = eps[i]->value(_chain);
		double const *mu = eps[i]->parents()[0]->value(_chain);
		for (unsigned int j = 0; j < m; ++j) {
		    eval[m*i + j] = mu[j] + (Y[j] - mu[j]) * sigma_ratio[j];
		}
	    }
	    _eps->setValue(eval, _chain);
	}
	
    }
}
