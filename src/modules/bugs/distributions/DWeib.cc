#include <config.h>
#include "DWeib.h"

#include <rng/RNG.h>
#include <util/nainf.h>

#include <vector>
#include <cmath>

#include <JRmath.h>

using std::vector;
using std::string;

#define SHAPE(par) (*par[0])
#define RATE(par) (*par[1])

namespace jags {
namespace bugs {

    DWeib::DWeib()
	: RScalarDist("dweib", 2, DIST_POSITIVE)
    {}
    
    string DWeib::alias() const
    {
	return "dweibull";
    }
    
    bool 
    DWeib::checkParameterValue (vector<double const *> const &par) const
    {
	return SHAPE(par) >= 0 && RATE(par) >  0;
    }
    
    double DWeib::d(double x, PDFType,
		    vector<double const *> const &par, bool give_log) const
	
    {
	return dweibull2(x, SHAPE(par), RATE(par), give_log);
    }

    double DWeib::p(double q, vector<double const *> const &par, bool lower, 
		    bool give_log) const
    {
	return pweibull2(q, SHAPE(par), RATE(par), lower, give_log);
    }

    double DWeib::q(double p, vector<double const *> const &par, bool lower, 
		    bool log_p) const
    {
	return qweibull2(p, SHAPE(par), RATE(par), lower, log_p);
    }

    double DWeib::r(vector<double const *> const &par, RNG *rng) const
    {
	return rweibull2(SHAPE(par), RATE(par), rng);
    }

}}
