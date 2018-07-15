#include <config.h>
#include <util/logical.h>
#include <util/nainf.h>
#include <util/dim.h>
#include <module/ModuleError.h>

#include "SumDist.h"

#include <cfloat>
#include <cmath>
#include <numeric>

using std::vector;
using std::fabs;
using std::sqrt;
using std::accumulate;

namespace jags {
namespace bugs {

    static const double TOL = sqrt(DBL_EPSILON);
    
    static double evaluate(vector <double const *> const &args,
			   vector<unsigned long> const &lengths)
    {
	double value = 0;
	for (unsigned long j = 0; j < args.size(); ++j) {
	    double const *a = args[j];
	    value = accumulate(a, a + lengths[j], value);
	}
	return value;
    }
    
    SumDist::SumDist()
	: VectorDist("sum", 0)
    {
    }
    
    bool SumDist::isDiscreteValued(vector<bool> const &mask) const
    {
	return allTrue(mask);
    }
    
    double SumDist::logDensity(double const *x, PDFType,
			       vector<double const *> const &par,
			       vector<unsigned long> const &lengths) const
    {
	return fabs(*x - evaluate(par, lengths)) > TOL ? JAGS_NEGINF : 0;
    }

    void SumDist::randomSample(double *x,
			       vector<double const *> const &par, 
			       vector<unsigned long> const &lengths,
			       RNG *) const
    {
	*x = evaluate(par, lengths);
    }

    bool SumDist::isSupportFixed(vector<bool> const &fixmask) const
    {
	return allTrue(fixmask);
    }

    unsigned long SumDist::df(vector<unsigned long> const &) const
    {
	return 0;
    }

    bool SumDist::checkParameterValue(vector<double const *> const &,
				      vector<unsigned long> const &) const
    {
	return true;
    }

    bool
    SumDist::checkParameterLength (vector<unsigned long> const &lengths) const
    {
	if (lengths.empty()) return false;
	for (unsigned long i = 1; i < lengths.size(); ++i) {
	    if (lengths[i] == 0)
		return false;
	}
	return true;
    }

    bool SumDist::checkParameterDiscrete(vector<bool> const &mask) const
    {
	for (unsigned long i = 1; i < mask.size(); ++i) {
	    if (mask[i] != mask[0])
		return false;
	}
	return true;
    }

    void SumDist::support(double *lower, double *upper,
			  vector<double const *> const &par,
			  vector<unsigned long> const &lengths) const
    {
	*lower = *upper = evaluate(par, lengths);
    }

    unsigned long SumDist::length(vector<unsigned long> const &) const
    {
	return 1;
    }

}}
