#include <config.h>
#include "DPick.h"

#include <util/nainf.h>
#include <rng/RNG.h>

#include <cmath>

using std::vector;
using std::min;
using std::max;
using std::log;

#define PROB(par) (*par[0])
#define VAL1(par) (*par[1])
#define VAL2(par) (*par[2])

namespace jags {
    namespace mix {

	DPick::DPick()
	    : ScalarDist("dpick", 3, DIST_SPECIAL)
	{}

	bool 
	DPick::checkParameterValue(vector<double const *> const &par) const
	{
	    if (PROB(par) <= 0.0 || PROB(par) >= 1.0) return false;

	    return true;
	}

	double DPick::logDensity(double x, 
			       PDFType type,
			       vector<double const *> const &par,
			       double const *lower, double const *upper) const
	{
	    double v1 = VAL1(par);
	    double v2 = VAL2(par);

	    if (v1 == x && v2 == x) {
		return 0;
	    }
	    else if (v1 == x) {
		return log(PROB(par));
	    }
	    else if (v2 == x) {
		return log(1 - PROB(par));
	    }
	    else {
		return JAGS_NEGINF;
	    }
	}

	double 
	DPick::randomSample(vector<double const *> const &par, 
			    double const *lower, double const *upper,
			    RNG *rng) const
	{
	    if (rng->uniform() <= PROB(par)) {
		return VAL1(par);
	    }
	    else {
		return VAL2(par);
	    }
	}

	double DPick::typicalValue(vector<double const *> const &par,
				   double const *lower,
				   double const *upper) const
	{
	    if (PROB(par) >= 0.5) {
		return VAL1(par);
	    }
	    else {
		return VAL2(par);
	    }
	}

	bool DPick::isSupportFixed(vector<bool> const &fixmask) const
	{
	    return fixmask[1] && fixmask[2];
	}

	void DPick::support(double *lower, double *upper, 
			    vector<double const *> const &par) const
	{
	    *lower = min(VAL1(par), VAL2(par));
	    *upper = max(VAL1(par), VAL2(par));
	}

	bool DPick::isDiscreteValued(std::vector<bool> const &mask) const
	{
	    return mask[1] && mask[2];
	}

    }
}
