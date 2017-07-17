#include <config.h>

#include "DOrdered.h"

#include <rng/RNG.h>
#include <util/nainf.h>

#include <cmath>

using std::log;
using std::string;
using std::vector;

#define MU(par) (*par[0])
#define CUT(par) (par[1])
#define NCUT(lengths) (lengths[1])

namespace jags {
    namespace glm  {

	DOrdered::DOrdered(string const &name)
	    : VectorDist(name, 2) 
	{}
	
	bool DOrdered::isDiscreteValued(vector<bool> const &mask) const
	{
	    return true;
	}
	
	bool
	DOrdered::checkParameterValue(vector<double const *> const &par,
				      vector<unsigned int> const &lengths) const
	{
	    double const * cut = CUT(par);
	    unsigned int ncut = NCUT(lengths);
	    
	    for (unsigned int i = 1; i < ncut; ++i) {
		if (cut[i] <= cut[i-1]) return false;
	    }
	    return true;
	}

	double DOrdered::density(double x, double mu, double const *cut,
				 int ncut, bool give_log) const
	{
	    int y = static_cast<int>(x) - 1;

	    if (y < 0 || y > ncut) {
		return JAGS_NEGINF;
	    }
	    else if (y == 0) {
		return p(cut[0], mu, true, give_log);
	    }
	    else if (y == ncut) {
		return p(cut[ncut - 1], mu, false, give_log);
	    }
	    else {
		double delta = p(cut[y], mu, true, false) -
		    p(cut[y-1], mu, true, false);
		return give_log ? log(delta) : delta;
	    }
	}
	
	double
	DOrdered::logDensity(double const *x, unsigned int length,
			     PDFType type,
			     vector<double const *> const &par,
			     vector<unsigned int> const &lengths,
			     double const *lower, double const *upper) const
	{
	    return density(*x, MU(par), CUT(par), NCUT(lengths), true);
	}

	void DOrdered::randomSample(double *x, unsigned int length,
				    vector<double const *> const &par,
				    vector<unsigned int> const &lengths,
				    double const *lower, double const *upper,
				    RNG *rng) const
	{
	    double y = r(MU(par), rng);
	    for (unsigned int i = 0; i < NCUT(lengths); ++i) {
		if (y <= CUT(par)[i]) {
		    *x = i + 1;
		    return;
		}
	    }
	    *x = NCUT(lengths) + 1;
	}

	void DOrdered::support(double *lower, double *upper,
			       unsigned int length,
			       vector<double const *> const &par,
			       vector<unsigned int> const &lengths) const
	{
	    *lower = 1;
	    *upper = NCUT(lengths) + 1;
	}
	
	void
	DOrdered::typicalValue(double *x, unsigned int length,
			       vector<double const *> const &par,
			       vector<unsigned int> const &lengths,
			       double const *lower, double const *upper) const
	{

	    double y = MU(par);
	    for (unsigned int i = 0; i < NCUT(lengths); ++i) {
		if (y <= CUT(par)[i]) {
		    *x = i + 1;
		    return;
		}
	    }
	    *x = NCUT(lengths);
	}
	
	bool DOrdered::isSupportFixed(vector<bool> const &fixmask) const
	{
	    return true;
	}
	
	bool DOrdered::checkParameterLength(vector<unsigned int> const &lengths)
	    const
	{
	    return lengths[0] == 1 && NCUT(lengths) > 0;
	}
	
	unsigned int DOrdered::length(vector<unsigned int> const &lengths) const
	{
	    return 1;
	}
	
	double DOrdered::KL(vector<double const *> const &par0,
			    vector<double const *> const &par1,
			    vector<unsigned int> const &lengths) const
	{
	    double psum0 = 0, psum1 = 0, y = 0;
	    for (unsigned int i = 0; i <= NCUT(lengths); ++i) {
		double p0 = density(i+1, MU(par0), CUT(par0), NCUT(lengths),
		    false);
		double p1 = density(i+1, MU(par1), CUT(par1), NCUT(lengths),
		    false);
		if (p0 == 0) {
		    psum1 += p1;
		}
		else if (p1 == 0) {
		    return JAGS_POSINF;
		}
		else {
		    y += p0 * (log(p0) - log(p1));
		    psum0 += p0;
		    psum1 += p1;
		}
	    }
	    y /= psum0;
	    y -= (log(psum0) - log(psum1));
	    return y;
	}

    }
}
