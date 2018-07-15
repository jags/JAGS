#include <config.h>
#include "DNT.h"

#include <cmath>

#include <JRmath.h>

using std::vector;

static inline double MU(vector<double const*> const &par) { return *par[0]; }
static inline double TAU(vector<double const*> const &par) { return *par[1]; }
static inline double DF(vector<double const*> const &par) { return *par[2]; }

static inline double SIGMA(vector<double const*> const &par) {
    //Standard deviation of numerator
    return 1/sqrt(TAU(par));
}
static inline double DELTA(vector<double const*> const &par) {
    //Non-centrality parameter of scaled x
    return MU(par)/SIGMA(par);
}

namespace jags {
    namespace bugs {

	DNT::DNT()
	    : RScalarDist("dnt", 3, DIST_UNBOUNDED)
	{}

	bool DNT::checkParameterValue (vector<double const *> const &par) const
	{
	    return (TAU(par) > 0 && DF(par) > 0 && fabs(DELTA(par)) <= 37.62);
	}

	double DNT::d(double x, PDFType,
		      vector<double const *> const &par, bool give_log) const
	{
	    x /= SIGMA(par);
	    if (give_log) {
		return dnt(x, DF(par), DELTA(par), 1) - log(SIGMA(par));
	    }
	    else {
		return dnt(x, DF(par), DELTA(par), 0) / SIGMA(par);
	    }
	}

	double DNT::p(double x, vector<double const *> const &par, bool lower, 
		      bool use_log) const
	{
	    return pnt(x/SIGMA(par), DF(par), DELTA(par), lower, use_log);
	}
	
	double DNT::q(double p, vector<double const *> const &par, bool lower, 
		      bool log_p) const
	{
	    return qnt(p, DF(par), DELTA(par), lower, log_p) * SIGMA(par);
	}
	
	double DNT::r(vector<double const *> const &par, RNG *rng) const
	{
	    double k = DF(par);
	    return rnorm(MU(par), SIGMA(par), rng)/sqrt(rchisq(k, rng)/k);
	}
    }
}
