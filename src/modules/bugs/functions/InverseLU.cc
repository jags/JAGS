#include <config.h>

#include <util/dim.h>
#include <util/integer.h>
#include <matrix/matrix.h>
#include <module/ModuleError.h>

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
	    bool can_invert = inverse_lu (value, args[0], dims[0][0]);
	    if (!can_invert) {
		throwFuncError(this, "Cannot invert matrix. It may be singular");
	    }
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

	bool InverseLU::hasGradient(unsigned long i) const
	{
	    return i == 0;
	}
	
	void InverseLU::gradient(double *grad, vector<double const *> const &args,
				 vector<vector<unsigned long>> const &dims,
				 unsigned long i) const
	{
	    unsigned long n = dims[0][0];
	    vector<double> y(n*n);
	
	    bool can_invert = inverse_lu (y.data(), args[0], n);
	    if (!can_invert) {
		throwFuncError(this, "Cannot invert matrix. It may be singular");
	    }

	    for (unsigned int i = 0; i < n; ++i) {
		for (unsigned int j = 0; j < n; ++j) {
		    for (unsigned int k = 0; k < n; ++k) {
			for (unsigned int l = 0; l < n; ++l) {
			    grad[i + n*(j + n*(k + n*l))] -= y[i + n*k] * y[l + n*j];
			}
		    }
		}
	    }

	}


    }
}
