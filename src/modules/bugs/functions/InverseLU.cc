#include <config.h>
#include "matrix.h"
#include <util/dim.h>
#include <util/integer.h>

#include "InverseLU.h"

using std::vector;

namespace jags {
    namespace bugs {

	InverseLU::InverseLU (): ArrayFunction ("inverse.lu", 1)
	{
	}

	void InverseLU::evaluate (double *value, vector<double const *> const &args,
				vector<vector<unsigned long> > const &dims) const
	{
	    inverse_lu (value, args[0], dims[0][0]);
	}

	bool 
	InverseLU::checkParameterDim (vector<vector<unsigned long> > const &dims) const
	{
	    return isSquareMatrix(dims[0]) || isScalar(dims[0]);
	}

	vector<unsigned long> 
	InverseLU::dim (vector<vector<unsigned long> > const &dims,
		      vector<double const *> const &values) const
	{
	    return dims[0];
	}

    }
}
