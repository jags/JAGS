/* Copyright (C) 2026 Marcel Jonker */

#include <config.h>

#include "EllipticalSlice.h"

#include <graph/StochasticNode.h>
#include <sampler/SingletonGraphView.h>
#include <rng/RNG.h>
#include <module/ModuleError.h>
#include <matrix/blas.h>
#include <matrix/matrix.h>
#include <util/integer.h>

#include <cmath>

using std::vector;
using std::cos;
using std::sin;
using std::isfinite;

namespace jags {
    namespace bugs {

	static const int one = 1;
	static const double TWO_PI = 6.283185307179586476925286766559;
	static const unsigned long MAX_SHRINK = 1000;

	/*
	 * Cholesky factorization L L^T of the precision matrix T. Returns
	 * false if T is not positive definite.
	 */
	static bool
	compute_chol(vector<double> &chol, double const *priorprec, unsigned long d)
	{
	    try {
		cholesky(chol.data(), priorprec, d);
	    }
	    catch (...) {
		return false;
	    }
	    return true;
	}

	/*
	 * Draw nu from N(0, Sigma), where Sigma is the inverse of the
	 * precision matrix T. If L L^T = T and z ~ N(0, I), then solving
	 * L^T nu = z gives nu covariance T^{-1}.
	 */
	static void
	draw_prior(vector<double> &nu,
		   vector<double> const &chol,
		   unsigned long d,
		   RNG *rng)
	{
	    for (unsigned long i = 0; i < d; ++i) {
		nu[i] = rng->normal();
	    }

	    int nr = asInteger(d);
	    jags_dtrsv("L", "T", "N", &nr, chol.data(), &nr, nu.data(), &one);
	}

	EllipticalSlice::EllipticalSlice(SingletonGraphView const *gv,
					 unsigned int chain)
	    : _gv(gv),
	      _chain(chain),
	      _length(gv->length()),
	      _xcur(_length),
	      _xprop(_length),
	      _nu(_length),
	      _chol(_length * _length),
	      _fixed_prec(gv->node()->parents()[1]->isFixed())
	{
	    gv->checkFinite(chain); //Check validity of initial values

	    if (_fixed_prec) {
		/* The precision matrix cannot change: factorize it once here instead of in every call to update. */
		double const *priorprec = gv->node()->parents()[1]->value(chain);
		if (!compute_chol(_chol, priorprec, _length)) {
		    throwNodeError(gv->node(), "Cannot calculate Cholesky decomposition of dmnorm precision matrix");
		}
	    }
	}

	void
	EllipticalSlice::update(RNG *rng)
	{
	    _gv->getValue(_xcur, _chain);

	    double ll_cur = _gv->logLikelihood(_chain);
	    if (!isfinite(ll_cur)) {
		if (ll_cur > 0) {
		    throwNodeError(_gv->node(), "Elliptical slice stuck at value with infinite likelihood");
		}
		else {
		    throwNodeError(_gv->node(), "Current value is inconsistent with data");
		}
	    }

	    /* Re-read parent values: they may be stochastic. */
	    StochasticNode const *snode = _gv->node();
	    double const *priormean = snode->parents()[0]->value(_chain);

	    if (!_fixed_prec) {
		double const *priorprec = snode->parents()[1]->value(_chain);
		if (!compute_chol(_chol, priorprec, _length)) {
		    throwNodeError(_gv->node(), "Cannot calculate Cholesky decomposition of dmnorm precision matrix");
		}
	    }

	    draw_prior(_nu, _chol, _length, rng);

	    /* Slice threshold on the log likelihood scale. */
	    double const log_y = ll_cur - rng->exponential();

	    /* Initial bracket contains theta = 0, which gives the current value. */
	    double theta = TWO_PI * rng->uniform();
	    double theta_min = theta - TWO_PI;
	    double theta_max = theta;

	    for (unsigned long s = 0; s < MAX_SHRINK; ++s) {
		double const ct = cos(theta);
		double const st = sin(theta);
		for (unsigned long i = 0; i < _length; ++i) {
		    _xprop[i] = priormean[i] + (_xcur[i] - priormean[i]) * ct + _nu[i] * st;
		}
		double ll_prop = 0.0;
		try {
		    _gv->setValue(_xprop, _chain);
		    ll_prop = _gv->logLikelihood(_chain);
		}
		catch (...) {
		    /* Restore the current value before rethrowing errors from descendants. */
		    _gv->setValue(_xcur, _chain);
		    throw;
		}
    
		if (isfinite(ll_prop) && ll_prop > log_y) {
		    return;
		}
    
		if (theta < 0.0) {
		    theta_min = theta;
		}
		else {
		    theta_max = theta;
		}
		theta = theta_min + (theta_max - theta_min) * rng->uniform();
	    }

	    /* Shrinkage budget exhausted. Restore the current value before reporting failure. */
	    _gv->setValue(_xcur, _chain);
	    throwNodeError(_gv->node(), "Elliptical slice shrinkage exhausted without acceptance");
	}

	bool
	EllipticalSlice::isAdaptive() const
	{
	    return false;
	}

	void
	EllipticalSlice::adaptOff()
	{
	}

	bool
	EllipticalSlice::checkAdaptation() const
	{
	    return true;
	}

    }
}
