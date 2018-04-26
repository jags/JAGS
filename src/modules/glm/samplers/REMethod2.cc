#include <config.h>

#include "REMethod2.h"
#include "GLMMethod.h"
#include "Outcome.h"

#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <JRmath.h>

#include <cmath>
#include <algorithm>
#include <set>

using std::vector;
using std::sqrt;
using std::fill;
using std::set;

extern cholmod_common *glm_wk;

namespace jags {
    namespace glm {

	static unsigned int sumLengths(vector<Outcome*> const &outcomes)
	{
	    unsigned int s = 0;
	    for (unsigned int i = 0; i < outcomes.size(); ++i) {
		s += outcomes[i]->length();
	    }
	    return s;
	}
	
	REMethod2::REMethod2(SingletonGraphView const *tau,
			     GLMMethod const *glmmethod)
	    : _tau(tau), _eps(glmmethod->_view),
	      _outcomes(glmmethod->_outcomes),
	      _x(glmmethod->_x), _chain(glmmethod->_chain)
	{
	    vector<StochasticNode*> const &enodes = _eps->nodes();
	    vector<StochasticNode*> const &schild = tau->stochasticChildren();

	    set<StochasticNode*> sset;
	    sset.insert(schild.begin(), schild.end());

	    for (unsigned int i = 0; i < enodes.size(); ++i) {
		if (sset.count(enodes[i])) {
		    if (tau->isDependent(enodes[i]->parents()[0])) {
			throwLogicError("Invalid REMethod2");
		    }
		    _indices.push_back(i);
		}
	    }

	    if (_indices.size() != schild.size()) {
		throwLogicError("Invalid REMethod2");
	    }

	    unsigned long nrow = sumLengths(_outcomes);
	    unsigned long ncol = tau->stochasticChildren()[0]->length();
	    _z = cholmod_allocate_dense(nrow, ncol, nrow, CHOLMOD_REAL,
					glm_wk);
	}

	REMethod2::~REMethod2()
	{
	    cholmod_free_dense(&_z, glm_wk);
	}
	
	void REMethod2::calDesignSigma()
	{
	    //Sanity checks
	    //unsigned int Neps = _eps->nodes().size();
	    if (_z->nrow != _x->nrow) {
		throwLogicError("Row mismatch in REMethod2");
	    }
	    /*
	    if (_x->ncol != _z->ncol * Neps || _x->ncol != _eps->length()) {
		throwLogicError("Column mismatch in REMethod2");
	    }
	    */
	    
	    //Set up access to sparse design matrix for eps
	    int const *Xp = static_cast<int const*>(_x->p);
	    int const *Xi = static_cast<int const*>(_x->i);
	    double const *Xx = static_cast<double const*>(_x->x);
	    
	    //Set up access to dense design matrix for sigma
	    double *Zx = static_cast<double*>(_z->x);

	    //Set all elements of _z to zero.
	    fill(Zx, Zx + _z->nzmax, 0);

	    //If there are m columns of _z then _z[,i] is the sum of
	    //every mth column of _x (starting with _x[,i]),
	    //multiplied by the corresponding random effect

	    vector<StochasticNode*> const &eps=_eps->nodes();
	    for (unsigned int j = 0; j < _indices.size(); ++j) {
		unsigned int i = _indices[j];
		double const *eval = eps[i]->value(_chain);
		double const *emean = eps[i]->parents()[0]->value(_chain);
		for (unsigned long zcol = 0; zcol < _z->ncol; ++zcol) {
		    unsigned long xcol = i * _z->ncol + zcol;
		    for (int xi = Xp[xcol]; xi < Xp[xcol+1]; ++xi) {
			int row = Xi[xi];
			int zi = _z->nrow * zcol + row;
			Zx[zi] += Xx[xi] * (eval[zcol] - emean[zcol]);
		    }
		}
	    }
	}

	void REMethod2::update(RNG *rng) {

	    updateTau(rng); //Sufficient parameterization
	    updateSigma(rng); //Ancillary parameterization
	}

	void REMethod2::calCoefSigma(double *A, double *b, double const *sigma0,
				     unsigned long m) const
	{
	    double const *Zx = static_cast<double const *>(_z->x);
	    unsigned long N = _outcomes.size();
	    
	    unsigned long xrow = 0;
	    for (unsigned long i = 0; i < N; ++i) {
		unsigned long n = _outcomes[i]->length();
		if (n == 1) {
		    //Scalar outcome
		    double Y = _outcomes[i]->value();
		    double mu = _outcomes[i]->mean();
		    double lambda = _outcomes[i]->precision();
		    vector<double> X(m);
		    for (unsigned int j = 0; j < m; ++j) {
			X[j] =  Zx[j * _z->nrow + xrow]/sigma0[j];
		    }
		    for (unsigned int j = 0; j < m; ++j) {
			for (unsigned int k = 0; k < m; ++k) {
			    A[j*m + k] += X[j] * X[k] * lambda;
			}
			b[j] += (Y - mu) * X[j] * lambda;
		    }
		}
		else {
		    //Multivariate outcome
		    double const *Y = _outcomes[i]->vvalue();
		    double const *mu = _outcomes[i]->vmean();
		    double const *lambda = _outcomes[i]->vprecision();

		    // delta = Y - mu
		    vector<double> delta(n);
		    for (unsigned int p = 0; p < n; ++p) {
			delta[p] = Y[p] - mu[p];
		    }
		    
		    //X = Scaled design matrix for outcome i
		    vector<double> X(m*n); //n x m matrix
		    for (unsigned long j = 0; j < m; ++j) {
			double const *Zj = &Zx[j * _z->nrow + xrow];
			for (unsigned int p = 0; p < n; ++p) {
			    X[j*n + p] = Zj[p]/sigma0[j];
			}
		    }

		    // Tx = lambda %*% X (n x m matrix)
		    vector<double> TX(n*m, 0); 
		    for (unsigned int j = 0; j < m; ++j) {
			for (unsigned int p = 0; p < n; ++p) {
			    b[j] += delta[p] * X[j*n + p];
			    for (unsigned int q = 0; q < n; ++q) {
				TX[j*n + p] += lambda[n*p + q] * X[j*n + q];
			    }
			}
		    }

		    for (unsigned int j = 0; j < m; ++j) {
			for (unsigned int p = 0; p < n; ++p) {
			    b[j] += delta[p] * TX[j*n + p];
			    for (unsigned int k = 0; k < m; ++k) {
				A[j*m + k] += X[j*n + p] * TX[k*n + p];
			    }
			}
		    }
		}
		xrow += n;
	    }
	}
	
	double REMethod2::logLikelihoodSigma(double const *sigma,
					    double const *sigma0,
					    unsigned int m) const
	{
	    vector<double> A(m*m, 0);
	    vector<double> b(m, 0);
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
	
    }
}
