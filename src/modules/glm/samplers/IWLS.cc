#include <config.h>
#include "IWLS.h"
#include <graph/LinkNode.h>
#include <graph/StochasticNode.h>
//#include <function/LinkFunction.h>
#include <distribution/Distribution.h>
#include <rng/RNG.h>

#include <stdexcept>
#include <cmath>

extern "C" {
#include <cholmod.h>
}

extern cholmod_common *glm_wk;

using std::string;
using std::vector;
using std::logic_error;
using std::log;
using std::exp;
using std::fabs;

static unsigned int nchildren(GraphView const *view)
{
    return view->stochasticChildren().size();
}

double logDet(cholmod_factor *F)
{
    if (!F->is_ll && !F->is_monotonic) {
	throw logic_error("Non-monotonic simplicial factor in logDet");
    }

    int *Fp = static_cast<int*>(F->p);
    double *Fx = static_cast<double*>(F->x);

    double y = 0;
    for (unsigned int r = 0; r < F->n; ++r) {
	y += log(Fx[Fp[r]]);
    }
    return F->is_ll ? 2*y : y;
}

#define MAX_ITER 100

namespace glm {


    IWLS::IWLS(GraphView const *view, 
	       vector<GraphView const *> const &sub_views,
	       unsigned int chain)
	: GLMMethod(view, sub_views, chain, true),
	  _link(nchildren(view)), _family(nchildren(view)), 
	  _init(true), _w(0)
    {
	vector<StochasticNode const*> const &children =
	    view->stochasticChildren();
	
	for (unsigned int i = 0; i < children.size(); ++i) {
	    _link[i] = dynamic_cast<LinkNode const*>(children[i]->parents()[0]);
	    _family[i] = getFamily(children[i]);
	}
    }
    
    string IWLS::name() const
    {
	return "IWLS";
    }

    double IWLS::getPrecision(unsigned int i) const
    {
	double w = _w;
	if (_family[i] == GLM_BINOMIAL) {
	    Node const *size = _view->stochasticChildren()[i]->parents()[1];
	    w *= size->value(_chain)[0];
	}

	double grad = _link[i]->grad(_chain);
	return (w * grad * grad)/ var(i);
    }

    double IWLS::getValue(unsigned int i) const
    {
	Node const *child = _view->stochasticChildren()[i];

	double y = child->value(_chain)[0];
	if (_family[i] == GLM_BINOMIAL) {
	    double N = child->parents()[1]->value(_chain)[0];
	    y /= N;
	}

	double mu = _link[i]->value(_chain)[0];
	double eta = _link[i]->eta(_chain);
	double grad = _link[i]->grad(_chain);

	return eta + (y - mu) / grad;
    }
 
    double IWLS::var(unsigned int i) const
    {
	double mu = _link[i]->value(_chain)[0];

	switch(_family[i]) {
	case GLM_BERNOULLI: case GLM_BINOMIAL:
	    return  mu * (1 - mu);
	    break;
	case GLM_POISSON:
	    return mu;
	    break;
	case GLM_NORMAL:
	    return 1;
	    break;
	case GLM_UNKNOWN:
	    throw logic_error("Unknown GLM family in IWLS");
	}
	return 0; //-Wall
    }

    double IWLS::logPTransition(vector<double> const &xorig, 
				vector<double> const &x,
				double *b, cholmod_sparse *A)
    {
	unsigned int n = _view->length();
	
	//Difference between new and old values
	cholmod_dense *delta = cholmod_allocate_dense(n, 1, n, CHOLMOD_REAL, 
						      glm_wk);
	double *dx = static_cast<double*>(delta->x);
	for (unsigned int i = 0; i < n; ++i) {
	    dx[i] = x[i] - xorig[i];
	}
	
	int ok = cholmod_factorize(A, _factor, glm_wk);
	if (!ok) {
	    throw logic_error("Cholesky decomposition failure in IWLS");
	}

	//Posterior mean
	cholmod_dense *mu = cholmod_solve(CHOLMOD_A, _factor, delta, glm_wk);
	double *mux = static_cast<double*>(mu->x);

	//Setup pointers to sparse matrix A
	int *Ap = static_cast<int*>(A->p);
	int *Ai = static_cast<int*>(A->i);
	double *Ax = static_cast<double*>(A->x);

	double deviance = 0;
	for (unsigned int r = 0; r < n; ++r) {
	    double Adr = 0;
	    for (int j = Ap[r]; j < Ap[r+1]; ++j) {
		Adr += Ax[j] * dx[Ai[j]];
	    }
	    deviance += dx[r] * (Adr - 2 * b[r])  + b[r] * mux[r];
	}
	deviance -= logDet(_factor);

	cholmod_free_dense(&delta, glm_wk);
	cholmod_free_dense(&mu, glm_wk);

	return -deviance/2;
    }
		
    void IWLS::update(RNG *rng)
    {
	if (_init) {
	    _w = 0;
	    for (unsigned int i = 0; i < MAX_ITER; ++i) {
		_w += 1.0/MAX_ITER;
		updateLM(rng, false);
            }
	    _init = false;
	}
	
	double *b1, *b2;
	cholmod_sparse *A1, *A2;
	double logp = 0;
	
	vector<double> xold(_view->length());
	_view->getValue(xold, _chain);
	calCoef(b1, A1);
	
	logp -= _view->logFullConditional(_chain);
	updateLM(rng);
	logp += _view->logFullConditional(_chain);

	vector<double> xnew(_view->length());
	_view->getValue(xnew, _chain);
	calCoef(b2, A2);

	logp -= logPTransition(xold, xnew, b1, A1);
	logp += logPTransition(xnew, xold, b2, A2);
	
	cholmod_free_sparse(&A1, glm_wk);
	cholmod_free_sparse(&A2, glm_wk);
	delete [] b1; delete [] b2;
	
	//Acceptance step
	if (rng->uniform() > exp(logp)) {
	    _view->setValue(xold, _chain);
	}
    }
}
