#include <config.h>

#include "HolmesHeldBlock.h"
#include "KS.h"

#include <graph/StochasticNode.h>
#include <graph/LinkNode.h>
#include <sampler/GraphView.h>
#include <rng/TruncatedNormal.h>
#include <rng/RNG.h>
#include <module/ModuleError.h>

extern "C" {
#include <cs.h>
}

#include <cmath>
#include <algorithm>

// Regularization penalty for precision 
#define REG_PENALTY 0.001

using std::vector;
using std::string;
using std::sqrt;
using std::copy;

extern cholmod_common *glm_wk;


/* FIXME: copied from AlbertChib */

//Random sample from a left-truncated logistic distribution
static double llogit(double left, RNG *rng, double mu)
{
    double qleft = 1/(1 + exp(mu-left));
    double x = qleft + (1 - qleft) * rng->uniform();
    return mu + log(x) - log(1 - x);
}

//Random sample from a right-truncated logistic distribution
static double rlogit(double right, RNG *rng, double mu)
{
    double qright = 1/(1 + exp(mu-right));
    double x = qright * rng->uniform();
    return mu + log(x) - log(1 - x);
}

namespace glm {
    
    HolmesHeldBlock::HolmesHeldBlock(GraphView const *view,
			   vector<GraphView const *> const &sub_views,
			   unsigned int chain)
	: BinaryGLM(view, sub_views, chain)
    {
    }
    
    void HolmesHeldBlock::updateAuxiliary(cholmod_dense *u, cholmod_factor *F, 
					  RNG *rng)
    {
	// Set up access to u
	double *ux = static_cast<double*>(u->x);

	// Set up access to _factor F
	int *fp = static_cast<int*>(F->p);
	double *fx = static_cast<double*>(F->x);

	//Transpose and permute the design matrix
	cholmod_sparse *t_x = cholmod_transpose(_x, 1, glm_wk);
        cholmod_sparse *Pt_x = cholmod_submatrix(t_x, F->Perm, F->n, NULL, -1, TRUE, TRUE, glm_wk);
        cholmod_free_sparse(&t_x, glm_wk);
	
	vector<StochasticNode const *> const &schildren =  
	    _view->stochasticChildren();
	unsigned int nrow = schildren.size();
	unsigned int ncol = _view->length();
	
	if (_x->nrow != nrow || _x->ncol != ncol) 
	    throwLogicError("Dimension mismatch in HolmesHeldBlock");

	// Setup vector of updates to RHS of equation L %*% u = b
	cholmod_dense *deltab = cholmod_zeros(ncol, 1, CHOLMOD_REAL, glm_wk);
	double *deltabx = static_cast<double*>(deltab->x);

	for (unsigned int r = 0; r < nrow; ++r) {

	    // FIXME: We could take a shallow copy here
	    int ir = r; 
	    cholmod_sparse *xr = cholmod_submatrix(Pt_x, 0, -1, &ir, 1, 
						   true, true, glm_wk);
	    double *xrx = static_cast<double *>(xr->x);
	    int const *xri = static_cast<int const*>(xr->i);

	    // Solve L %*% D %*% v = xr
	    cholmod_sparse *v = cholmod_spsolve(CHOLMOD_LD, F, xr, glm_wk);
	    int const *vi = static_cast<int const *>(v->i);
	    double const *vx = static_cast<double const *>(v->x);

	    // Calculate mean and variance of linear predictor
	    double eta_mean = getMean(r); 
	    double eta_var = 0;
	    for (unsigned int j = 0; j < v->nrow; ++j) {
		eta_mean += vx[j] * ux[vi[j]];
		if (F->is_ll) {
		    // LL' decomposition
		    eta_var += vx[j] * vx[j];
		}
		else {
		    // LDL' decomposition. The diagonal D matrix is
		    // stored in the first element of each column of F
		    eta_var += vx[j] * vx[j] * fx[fp[vi[j]]];
		}
	    }
	    cholmod_free_sparse(&v, glm_wk);

	    // Sample linear predictor
	    // FIXME: It's not this simple! Since we are sampling
	    // z[r] at the same time, we need to account for the fact that
	    // z[r] has an interval restriction (<=0 or >=0) when sampling
	    // eta.  This requires something like a rejection sampling step
	    double eta = eta_mean + sqrt(eta_var) * rng->normal();

	    // Now sample auxiliary variable from truncated
	    // distribution given eta
	    double yr = *schildren[r]->value(_chain);
	    switch(_outcome[r]) {
	    case BGLM_NORMAL:
		_z[r] = yr;
		break;
	    case BGLM_PROBIT:
		if (yr == 1) {
		    _z[r] = lnormal(0, rng, eta_mean, sqrt(eta_var + 1));
		}
		else if (yr == 0) {
		    _z[r] = rnormal(0, rng, eta_mean, sqrt(eta_var + 1));
		}
		else {
		    throwLogicError("Invalid child value in HolmesHeldBlock");
		}
		break;
	    case BGLM_LOGIT: 
		if (yr == 1) {
		    _z[r] = llogit(0, rng, eta); //FIXME
		}
		else if (yr == 0) {
		    _z[r] = rlogit(0, rng, eta); //FIXME
		}
		else {
		    throwLogicError("Invalid child value in HolmesHeldBlock");
		}
                _tau[r] = REG_PENALTY + 1/sample_lambda(_z[r] - eta, rng);
		break;
	    }

	    // Prepare for update
	    double tau_r = getPrecision(r);
	    double sigma_r = sqrt(tau_r);
	    double delta_r = tau_r * (_z[r] - getMean(r));
	    
	    for (unsigned int j = 0; j < xr->nrow; ++j) {
		deltabx[xri[j]] = xrx[j] * delta_r;
		xrx[j] *= sigma_r;
	    }
	    
	    // Update factorization
	    // deltab should be zero on exit
	    cholmod_updown_solve(true, xr, F, u, deltab, glm_wk);

	    // Debuggin
	    for (unsigned int i = 0; i < ncol; ++i) {
		if (deltabx[i] != 0) throwLogicError("deltab not reset to zero");
	    }

	    cholmod_free_sparse(&xr, glm_wk);
	}

	cholmod_free_dense(&deltab, glm_wk);
	cholmod_free_sparse(&Pt_x, glm_wk);
    }

    string HolmesHeldBlock::name() const
    {
	return "Holmes-Held Block";
    }

    void HolmesHeldBlock::update(RNG *rng) 
    {
	updateLM(rng, 0);
    }

}
    
