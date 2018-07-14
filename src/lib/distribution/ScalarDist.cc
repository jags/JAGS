#include <config.h>
#include <distribution/ScalarDist.h>
#include <util/nainf.h>
#include <util/dim.h>

#include <stdexcept>
#include <cmath>
#include <algorithm>

using std::max;
using std::min;
using std::string;
using std::vector;
using std::length_error;
using std::logic_error;
using std::count_if;

namespace jags {

ScalarDist::ScalarDist(string const &name, unsigned int npar, Support support)
  : Distribution(name, npar), _support(support)
{
}

double ScalarDist::l(vector<double const *> const &) const
{
    double lb = JAGS_POSINF;
    
    switch(_support) {
    case DIST_UNBOUNDED:
	lb = JAGS_NEGINF;
	break;
    case DIST_POSITIVE: case DIST_PROPORTION:
	lb = 0;
	break;
    case DIST_SPECIAL:
	//You must overload this function 
	throw logic_error("Cannot call ScalarDist::l for special distribution");
    }
    
    return lb;
}

double ScalarDist::u(vector<double const *> const &) const
{
    double ub = JAGS_NEGINF;
    
    switch(_support) {
    case DIST_UNBOUNDED: case DIST_POSITIVE:
	ub = JAGS_POSINF;
	break;
    case DIST_PROPORTION:
	ub = 1.0;
	break;
    case DIST_SPECIAL:
	//You must overload this function 
	throw logic_error("Cannot call ScalarDist::u for special distribution");
    }
    return ub;
}

bool ScalarDist::isSupportFixed(vector<bool> const &) const
{
    if (_support == DIST_SPECIAL) {
	//You must overload this function 
	throw logic_error("Cannot call ScalarDist::isSupportFixed for special distribution");
    }
    else {
	return true;
    }
}

unsigned long ScalarDist::df() const
{
    return 1;
}

    double ScalarDist::KL(vector<double const *> const &par1,
			  vector<double const *> const &par2,
			  double const *lower, double const *upper,
			  RNG *rng, unsigned int nrep) const
    {
	double div = 0;
	
	for (unsigned int r = 0; r < nrep; ++r) {
	    double v1 = randomSample(par1, lower, upper, rng);
	    div += logDensity(v1, PDF_FULL, par1, lower, upper);
	    div -= logDensity(v1, PDF_FULL, par2, lower, upper);
	}
	return div / nrep;
    }

    double ScalarDist::KL(vector<double const *> const &,
			  vector<double const *> const &) const
    {
	return JAGS_NA;
    }
    
} //namespace jags
