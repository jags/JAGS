#include <config.h>

#include "REGamma.h"
#include "Outcome.h"

#include <JRmath.h>
#include <sampler/SingletonGraphView.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>

#include <cmath>
#include <algorithm>

using std::vector;
using std::sqrt;
using std::fill;

extern cholmod_common *glm_wk;

namespace jags {
    namespace glm {

	REGamma::REGamma(SingletonGraphView const *tau,
			 GraphView const *eps,
			 vector<SingletonGraphView const *> const &sub_eps,
			 vector<Outcome *> const &outcomes,
			 unsigned int chain)
	    : REMethod(tau, eps, sub_eps, outcomes, chain)
	{
	    unsigned int nrow = _outcomes.size();
	    unsigned int ncol = eps->nodes()[0]->length();
	    _z = cholmod_allocate_dense(nrow, ncol, nrow, CHOLMOD_REAL,
					glm_wk);

	    /* FIXME: if _z is allocated in the parent class then we don't
	       need to dynamically allocate the slicer */
	    vector<Node const *> const &par = tau->node()->parents();
	    double const *shape = par[0]->value(chain);
	    double const *rate = par[1]->value(chain);
	    double sigma = 1/sqrt(*tau->node()->value(chain));
	    _slicer = new REGammaSlicer(outcomes, _x, _z, shape, rate, sigma);
	}

	REGamma::~REGamma()
	{
	    cholmod_free_dense(&_z, glm_wk);
	    delete _slicer;
	}

	void REGamma::updateTau(RNG *rng)
	{
	    vector<Node const*> const &par = _tau->node()->parents();
	    double shape = *par[0]->value(_chain); 
	    double rate = *par[1]->value(_chain); //(1/scale)
    
	    // Likelihood
	    //vector<StochasticNode *> const &sch = _tau->stochasticChildren();
	    vector<StochasticNode *> const &eps = _eps->nodes();
	    for (unsigned int i = 0; i < eps.size(); ++i) {
		double Y = *eps[i]->value(_chain);
		shape += 0.5;
		rate += Y * Y / 2.0;
	    }

	    double x = rgamma(shape, 1.0/rate, rng);
	    _tau->setValue(&x, 1, _chain);  
	}

	void REGamma::calDesignSigma()
	{
	    //Sanity checks
	    unsigned int Neps = _eps->nodes().size();
	    if (_x->nrow != _outcomes.size() || _z->nrow != _x->nrow) {
		throwLogicError("Row mismatch in REGamma");
	    }
	    if (_x->ncol != _z->ncol * Neps || _x->ncol != _eps->length()) {
		throwLogicError("Column mismatch in REGamma");
	    }
	    
	    //Get current values of random effects
	    vector<double> eps(_eps->length());
	    _eps->getValue(eps, _chain);
	    
	    //Set up access to sparse design matrix for eps
	    int *Xp = static_cast<int*>(_x->p);
	    int *Xi = static_cast<int*>(_x->i);
	    double *Xx = static_cast<double*>(_x->x);
	    
	    //Set up access to dense design matrix for sigma
	    double *Zx = static_cast<double*>(_z->x);

	    //Set all elements of _z to zero.
	    fill(Zx, Zx + _z->nzmax, 0);

	    /*
	      std::cout << "_z->ncol = " << _z->ncol << ", _z->nrow = "
	      << _z->nrow << std::endl;
	      std::cout << "_x->ncol = " << _x->ncol << ", _x->nrow = "
	      << _x->nrow << std::endl;
	      std::cout << "Neps = " << Neps << std::endl;
	    */
	    
	    /*
	      for (unsigned int zc = 0; zc < _z->ncol; ++zc) {
	      for (unsigned int i = 0; i < Neps; ++i) {
	      int xc = i * _z->ncol + zc;
	      for (int r = Xp[xc]; r < Xp[xc+1]; ++r) {
	      int xi = Xi[r];
	      int zi = _z->nrow * zc + r;
	      //_z[r,zc] += _x[r,xc] * exp[xc]
	      std::cout << "zc = " << zc << ", i = " << i
	      << ", xc = " << xc << ", xi = " << xi
	      << ", zi = "  << zi
	      << ", r = " << r << std::endl;
	      Zx[zi] += Xx[xi] * eps[xc];
	      }
	      }
	      }
	    */

	    //If there are m columns of _z then _z[,i] is the sum of
	    //every mth column of _x (starting with _x[,i]),
	    //multiplied by the corresponding random effect
	    for (unsigned int zcol = 0; zcol < _z->ncol; ++zcol) {
		for (unsigned int i = 0; i < Neps; ++i) {
		    int xcol = i * _z->ncol + zcol;
		    for (int xi = Xp[xcol]; xi < Xp[xcol+1]; ++xi) {
			int row = Xi[xi];
			int zi = _z->nrow * zcol + row;
			//_z[row,zcol] += _x[row,xcol] * exp[xcol]
			Zx[zi] += Xx[xi] * eps[xcol];
		    }
		}
	    }

	}

	void REGamma::updateSigma(RNG *rng)
	{
	    double tau = _tau->node()->value(_chain)[0];
	    double sigma0 = 1/sqrt(tau);

	    calDesignSigma();
	    
	    _slicer->setSigma(sigma0);
	    _slicer->update(rng);
	    double sigma1 = _slicer->value();

	    //Set new value of precision parameter
	    double x = 1/(sigma1 * sigma1);
	    _tau->setValue(&x, 1, _chain);

	    //Rescale random effects
	    vector<double> eps(_eps->length());
	    _eps->getValue(eps, _chain);
	    double sigma_ratio = sigma1/sigma0;
	    for (unsigned int i = 0; i < eps.size(); ++i) {
		eps[i] *= sigma_ratio;
	    }
	    _eps->setValue(eps, _chain);
	}

	void REGamma::update(RNG *rng) {
	    
	    // Update outcomes
	    for (vector<Outcome*>::const_iterator p = _outcomes.begin();
		 p != _outcomes.end(); ++p)
	    {
		(*p)->update(rng);
	    }
	    
	    updateEps(rng);
	    updateTau(rng); //Sufficient parameterization
	    updateSigma(rng); //Ancillary parameterization
	    updateTau(rng); //Sufficient parameterization
	}

	bool REGamma::isAdaptive() const
	{
	    return true;
	}
	
	void REGamma::adaptOff()
	{
	    _slicer->adaptOff();
	}
	
	bool REGamma::checkAdaptation() const
	{
	    return _slicer->checkAdaptation();
	}


    }
}
