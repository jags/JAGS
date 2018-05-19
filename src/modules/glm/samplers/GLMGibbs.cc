#include <config.h>

#include "GLMGibbs.h"
#include "Outcome.h"

#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <rng/TruncatedNormal.h>
#include <rng/RNG.h>
#include <module/ModuleError.h>

#include <cmath>

using std::vector;
using std::sqrt;

extern cholmod_common *glm_wk;

namespace jags {

namespace glm {

    GLMGibbs::GLMGibbs(GraphView const *view, 
			 vector<SingletonGraphView const *> const &sub_views,
			 vector<Outcome *> const &outcomes,
			 unsigned int chain)
	: GLMMethod(view, sub_views, outcomes, chain)
    {
	if (_view->length() != _sub_views.size()) {
	    throwLogicError("updateLMGibbs can only act on scalar nodes");
	}
	calDesign();
    }

    void GLMGibbs::update(RNG *rng) 
    {
	// Update element-wise. Less efficient than GLMBlock#update
	// but does not require a Cholesky decomposition, and is
	// necessary for truncated parameters

	// Update outcomes
	for (vector<Outcome*>::const_iterator p = _outcomes.begin();
	     p != _outcomes.end(); ++p)
	{
	    (*p)->update(rng);
	}
	
	double *b = 0;
	cholmod_sparse *A = 0;
	calCoef(b, A);

	if (A->stype != 0) {
	    throwLogicError("Wrong stype in GLMGibbs::update");
	}
	
	int nrow = _view->length();
	vector<double> theta(nrow);
	_view->getValue(theta, _chain);

	int *Ap = static_cast<int*>(A->p);
	int *Ai = static_cast<int*>(A->i);
	double *Ax = static_cast<double*>(A->x);

	//Extract diagonal from A
	vector<double> diagA(nrow);
	for (int c = 0; c < nrow; ++c) {
	    for (int j = Ap[c]; j < Ap[c+1]; ++j) {
		if (Ai[j] == c) {
		    diagA[c] = Ax[j];
		    break;
		}
	    }
	}

	//Update element-wise
	for (int i = 0; i < nrow; ++i) {
	    
	    double theta_old = theta[i];
		
	    double mu  = theta[i] + b[i]/diagA[i];
	    double sigma = sqrt(1/diagA[i]);
	    StochasticNode const *snode = _sub_views[i]->node();
	    double const *l = snode->lowerLimit(_chain);
	    double const *u = snode->upperLimit(_chain);
		
	    if (l && u) {
		theta[i] = inormal(*l, *u, rng, mu, sigma);
	    }
	    else if (l) {
		theta[i] = lnormal(*l, rng, mu, sigma);
	    }
	    else if (u) {
		theta[i] = rnormal(*u, rng, mu, sigma);
	    }
	    else {
		theta[i] = mu + rng->normal() * sigma;
	    }
		
	    double delta = theta[i] - theta_old;
	    for (int j = Ap[i]; j < Ap[i+1]; ++j) {
		b[Ai[j]] -= delta * Ax[j];
	    }
	}

	#pragma omp critical 
	cholmod_free_sparse(&A, glm_wk);
	delete [] b;
	
	_view->setValue(theta,  _chain);
    }

}}
