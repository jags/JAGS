#include <config.h>
#include "DNChisqr.h"

#include <JRmath.h>

using std::vector;
using std::string;

#define DF(par) (*par[0])
#define NCP(par) (*par[1])

namespace jags {
namespace bugs {

    DNChisqr::DNChisqr()
	: RScalarDist("dnchisqr", 2, DIST_POSITIVE)
    {}

    string DNChisqr::alias() const
    {
	return "dnchisq";
    }

    bool 
    DNChisqr::checkParameterValue (vector<double const *> const &par) const
    {
	return (DF(par) > 0 && NCP(par) >= 0);
    }

    double 
    DNChisqr::d(double x, PDFType,
		vector<double const *> const &par, bool give_log) const
    {
	return dnchisq(x, DF(par), NCP(par), give_log);
    }

    double 
    DNChisqr::p(double q, vector<double const *> const &par, bool lower, bool log_p)
	const
    {
	return pnchisq(q, DF(par), NCP(par), lower, log_p);
    }

    double DNChisqr::q(double p, vector<double const *> const &par, 
		       bool lower, bool log_p) const
    {
	return qnchisq(p, DF(par), NCP(par), lower, log_p);
    }

    double DNChisqr::r(vector<double const *> const &par, RNG *rng) const
    {
	return rnchisq(DF(par), NCP(par), rng);
    }

}}
