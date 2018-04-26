#include <config.h>
#include "InterpLin.h"

#include <util/nainf.h>
#include <util/dim.h>

using std::vector;

namespace jags {
namespace bugs {

    InterpLin::InterpLin() : ScalarVectorFunction("interp.lin", 3)
    {}
    
    double InterpLin::scalarEval(vector<double const *> const &args,
				 vector<unsigned long> const &lengths) const
    {
	double const *x = args[1];
	double const *y = args[2];
	double xnew = args[0][0];
	unsigned long N = lengths[1];
    
	if (xnew < x[0])
	    return y[0];
	else if (xnew >= x[N-1])
	    return y[N-1];
	else {
	    unsigned long i = 0;
	    for (; i < N-1; ++i) {
		if (xnew >= x[i] && xnew < x[i+1]) {
		    break;
		}
	    }
	    if (i == N-1) {
		/* This shouldn't happen, but we must guard against an
		   attempt to access an illegal element of x or y */
		return JAGS_NAN;
	    }
	    else {
		return y[i] + 
		    (xnew - x[i]) * (y[i+1] - y[i]) / (x[i+1] - x[i]);
	    }
	}
    }

    bool InterpLin::checkParameterLength(vector<unsigned long> const &lengths)
	const
    {
	return lengths[0] == 1 && lengths[1] != 0 && lengths[2] == lengths[1];
    }

    bool 
    InterpLin::checkParameterValue(vector <double const *> const &args,
				   vector <unsigned long> const &lengths) const
    {
        unsigned long N = lengths[1];
	for (unsigned long i = 1; i < N; ++i) {
	    if (args[1][i] <= args[1][i-1])
		return false;
	}
	return true;
    }

}}
