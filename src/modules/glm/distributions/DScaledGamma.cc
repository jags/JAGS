#include <config.h>
#include "DScaledGamma.h"

#include <JRmath.h>
#include <util/nainf.h>

using std::vector;

static inline double S(vector<double const *> const &par) {
    return (*par[0]) * (*par[0]);
}
static inline double logS(vector<double const *> const &par) {
    return 2 * log(*par[0]);
}
static inline double DF(vector<double const *> const &par) {
    return *par[1];
}

namespace jags {
    namespace glm {
	
	DScaledGamma::DScaledGamma()
	    : RScalarDist("dscaled.gamma", 2, DIST_POSITIVE)
	{}
	
	bool
	DScaledGamma::checkParameterValue (vector<double const *> const &par)
	    const
	{
	    return (*par[0] > 0 && *par[1] > 0);
	}

	double DScaledGamma::d(double x, PDFType type,
			       vector<double const *> const &par,
			       bool give_log) const
	{
	    if (give_log) {
		return dF(S(par) * x, DF(par), 1, true) + logS(par);
	    }
	    else {
		return dF(S(par) * x, DF(par), 1, false) * S(par);
	    }
	}

	double DScaledGamma::p(double x, vector<double const *> const &par,
			       bool lower, bool use_log) const
	{
	    return pF(S(par) * x, DF(par), 1, lower, use_log);
	}
	
	double DScaledGamma::q(double p, vector<double const *> const &par,
			       bool lower, bool log_p) const
	{
	    return qF(p, DF(par), 1, lower, log_p) / S(par);
	}
	
	double DScaledGamma::r(vector<double const *> const &par,
			       RNG *rng) const
	{
	    return rF(DF(par), 1, rng) / S(par);
	}
	
    }
}
