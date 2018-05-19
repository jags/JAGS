#include <config.h>

#include "HolmesHeldGibbs.h"
#include "Outcome.h"

#include <graph/StochasticNode.h>
#include <graph/LinkNode.h>
#include <sampler/GraphView.h>
#include <rng/TruncatedNormal.h>
#include <rng/RNG.h>
#include <module/ModuleError.h>
#include <sampler/SingletonGraphView.h>

#include <cmath>

using std::vector;
using std::sqrt;

extern cholmod_common *glm_wk;

namespace jags {
    namespace glm {
	
	HolmesHeldGibbs::HolmesHeldGibbs(GraphView const *view,
		 vector<SingletonGraphView const *> const &sv,
		 vector<Outcome *> const &outcomes,
		 unsigned int chain)
	    : GLMMethod(view, sv, outcomes, chain)
	{
	}

	void HolmesHeldGibbs::update(RNG *rng) 
	{
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
		throwLogicError("Wrong stype in HolmesHeldGibbs::update");
	    }
	    
	    int nrow = _view->length();
	    vector<double> theta(nrow);
	    _view->getValue(theta, _chain);

	    vector<double> theta_old(theta);
	    
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

	    //Transpose design matrix
	    cholmod_sparse *t_x = 0;
	    #pragma omp critical
	    t_x = cholmod_transpose(_x, 1, glm_wk);
	
	    double *xx = static_cast<double*>(t_x->x);
	    int *xp = static_cast<int*>(t_x->p);
	    int *xi = static_cast<int*>(t_x->i);

	    //Update element-wise
	    for (unsigned int c = 0; c <  t_x->ncol; ++c) {

		double eta = _outcomes[c]->mean();
		double z_old = _outcomes[c]->value();
		
		for (int r = xp[c]; r < xp[c+1]; ++r) {

		    int i = xi[r];
		    
		    double theta_i = theta[i];
		    
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

		    double delta = theta[i] - theta_i;
		    for (int j = Ap[i]; j < Ap[i+1]; ++j) {
			b[Ai[j]] -= delta * Ax[j];
		    }

		    eta += (theta[i] - theta_old[i]) * xx[r];
		}

		_outcomes[c]->update(eta, 0, rng);
		double z_new = _outcomes[c]->value();
		double tau_c = _outcomes[c]->precision();
		double delta_z = (z_new - z_old) * tau_c;
		
		for (int j = xp[c]; j < xp[c+1]; ++j) {
		    b[xi[j]] += delta_z * xx[j];
		}
	    }

	    #pragma omp critical
	    cholmod_free_sparse(&A, glm_wk);
	    delete [] b;
	    
	    _view->setValue(theta,  _chain);
	}

    }
}
    
