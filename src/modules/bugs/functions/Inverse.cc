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
			   vector<vector<unsigned long> > const &dims,
			   unsigned long i) const
    {
	unsigned long N = dims[0][0] * dims[0][1];
	vector<double> y(N);
	
	bool can_invert = inverse_chol (y.data(), args[0], dims[0][0]);
	if (!can_invert) {
	    throwFuncError(this, "Cannot invert matrix. It may not be positive definite");
	}
	
	for (unsigned int i = 0; i < N; ++i) {
	    for (unsigned int j = 0; j < N; ++j) {
		grad[i + j*N] = y[i] * y[j];
	    }
	}
    }


}}
