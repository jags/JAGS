#include <config.h>
#include <util/logical.h>
#include <util/integer.h>

#include "InProd.h"

#include "lapack.h"

using std::vector;

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
