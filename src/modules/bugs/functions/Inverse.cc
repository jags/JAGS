#include <config.h>

#include <util/dim.h>
#include <util/integer.h>
#include <matrix/matrix.h>
#include <module/ModuleError.h>

#include "Inverse.h"

using std::vector;
using std::string;

namespace jags {
namespace bugs {

    Inverse::Inverse (): ArrayFunction ("inverse.chol", 1)
    {
    }

    string Inverse::alias() const
    {
        return "inverse";
    }

    void Inverse::evaluate (double *value, vector<double const *> const &args,
			    vector<vector<unsigned long> > const &dims) const
    {
	bool can_invert = inverse_chol (value, args[0], dims[0][0]);
	if (!can_invert) {
	    throwFuncError(this, "Cannot invert matrix. It may not be positive definite");
	}
    }

    bool 
    Inverse::checkParameterDim (vector<vector<unsigned long> > const &dims) const
    {
	return isSquareMatrix(dims[0]) || isScalar(dims[0]);
    }

    vector<unsigned long> 
    Inverse::dim (vector<vector<unsigned long> > const &dims,
		  vector<double const *> const &) const
    {
	return dims[0];
    }

    bool Inverse::hasGradient(unsigned long i) const
    {
	return i == 0;
    }

    void Inverse::gradient(double *grad, vector<double const *> const &args,
			   vector<vector<unsigned long>> const &dims,
			   unsigned long i) const
    {
	unsigned long n = dims[0][0];
	vector<double> y(n*n);
	
	bool can_invert = inverse_chol (y.data(), args[0], n);
	if (!can_invert) {
	    throwFuncError(this, "Cannot invert matrix. It may not be positive definite");
	}

	/* The gradient for the Cholesky inverse is complicated by two factors
	   1) JAGS only uses the lower triangle of the input matrix. Hence the gradient with
	   respect to the upper triangle (l > k) is zero.
	   2) In order to preserve symmetry, the perturbation being tested is symmetric for
	   the off diagonal elements.
	*/
	for (unsigned int i = 0; i < n; ++i) {
	    for (unsigned int j = 0; j <= i; ++j) {
		for (unsigned int k = 0; k < n; ++k) {
		    for (unsigned int l = 0; l <= k; ++l) {
			double delta = y[i + n*k] * y[j + n*l] + y[i + n*l] * y[j + n*k];
			if (i == j) {
			    delta /= 2;
			}
			if (l == k) {
			    delta /= 2;
			}
			grad[i + n*(j + n*(k + n*l))] -= delta;
			grad[j + n*(i + n*(k + n*l))] -= delta;
		    }
		}
	    }
	}

    }

}}
