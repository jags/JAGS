#include <config.h>
#include "DGamPois.h"
#include <util/nainf.h>

#include <cmath>
#include <vector>
#include <algorithm>

#include <JRmath.h>

using std::vector;
using std::max;
using std::string;

#define MEAN(par) (*par[0])
#define SHAPE(par) (*par[1])

namespace jags {
namespace bugs {

DGamPois::DGamPois()
    : RScalarDist("dgampois", 2, DIST_POSITIVE, true)
{}

bool 
DGamPois::checkParameterValue (vector<double const *> const &par) const
{
	double m = MEAN(par);
	double s = SHAPE(par);
    double prob = s / (s + m);
    return (s >= 0 && prob > 0 && prob <= 1);
}

double
DGamPois::d(double x, PDFType type,
	   vector<double const *> const &par, bool give_log) 
    const
{
    if (SHAPE(par) == 0) {
	if (give_log) {
	    return (x == 0) ? 0 : JAGS_NEGINF;
	}
	else {
	    return (x == 0) ? 1 : 0;
	}
    }
    else {
		double m = MEAN(par);
		double s = SHAPE(par);
	    double prob = s / (s + m);
		return dnbinom(x, s, prob, give_log);
    }
}

double
DGamPois::p(double q, vector<double const *> const &par, bool lower, 
	   bool give_log) const
{
    if (SHAPE(par) == 0) {
	return give_log ? 0 : 1;
    }
    else {
		double m = MEAN(par);
		double s = SHAPE(par);
	    double prob = s / (s + m);
		return pnbinom(q, s, prob, lower, give_log);
    }
}

double 
DGamPois::q(double p, vector<double const *> const &par, bool lower, 
	   bool log_p) const
{
    if (SHAPE(par) == 0) {
	return 0;
    }
    else {
		double m = MEAN(par);
		double s = SHAPE(par);
	    double prob = s / (s + m);
		return qnbinom(p, s, prob, lower, log_p);
    }
}

double DGamPois::r(vector<double const *> const &par, RNG *rng) const
{
    if (SHAPE(par) == 0) {
	return 0;
    }
    else {
		double m = MEAN(par);
		double s = SHAPE(par);
	    double prob = s / (s + m);
		return rnbinom(s, prob, rng);
    }
}

double DGamPois::KL(vector<double const *> const &par0,
		   vector<double const *> const &par1) const
{
	double m0 = MEAN(par0);
	double s0 = SHAPE(par0);
    double p0 = s0 / (s0 + m0);
	double m1 = MEAN(par1);
	double s1 = SHAPE(par1);
    double p1 = s1 / (s1 + m1);

    if (fabs(s0 - s1) > 1e-16) {
	//We can't calculat Kullback-Leibler divergence in closed form when
	//s0 and s1 are different
	return JAGS_NA;
    }
    
    return s0 * (log(p0) - log(p1))  + 
	(1 - p0) * s0 * (log(1 - p0) - log(1 - p1)) / p0;
}

}}
