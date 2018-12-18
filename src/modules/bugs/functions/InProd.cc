#include <config.h>
#include <util/logical.h>
#include <util/integer.h>

#include "InProd.h"
#include "lapack.h"

#include <algorithm>

using std::vector;
using std::copy;

namespace jags {
namespace bugs {

    InProd::InProd () : ScalarVectorFunction ("inprod", 2)
    {
    }

    double InProd::scalarEval(vector<double const *> const &args,
			      vector<unsigned long> const &lengths) const
    {
        int one = 1, N = asInteger(lengths[0]);
        return F77_DDOT(&N, args[0], &one, args[1], &one);
    }

    bool InProd::isDifferentiable(unsigned long i) const
    {
	return i < 2;
    }
    
    void InProd::gradient(double *grad,
			  vector <double const *> const &args,
			  vector<unsigned long> const &lengths,
			  unsigned long i) const
    {
	unsigned long k = (i == 0) ? 1 : 0;

	for (unsigned long j = 0; j < lengths[k]; ++j) {
	    grad[j] += args[k][j];
	}
    }

    bool 
    InProd::checkParameterLength (vector<unsigned long> const &lengths) const
    {
	return (lengths[0] > 0) && (lengths[0] == lengths[1]);
    }

    bool InProd::isDiscreteValued(vector<bool> const &mask) const
    {
	return allTrue(mask);
    }

    bool 
    InProd::isScale(vector<bool> const &mask, vector<bool> const &fix) const
    {
	//Test for quadratic term
	if (mask[0] && mask[1])
	    return false;

	if (fix.empty()) {
	    return true;
	}
        else {
	    return (mask[0] || fix[0]) && (mask[1] || fix[1]); 
        }
    }

}}
