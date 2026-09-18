#include <config.h>
#include <util/dim.h>

#include "Transpose.h"

using std::vector;

static inline unsigned long NROW(vector<unsigned long> const &dim)
{
    return dim[0];
}

static inline unsigned long NCOL(vector<unsigned long> const &dim)
{
    return dim.size() > 1 ? dim[1] : 1;
}


namespace jags {
namespace bugs {

    Transpose::Transpose()
	: ArrayFunction("t",1)
    {
    }

    void
    Transpose::evaluate (double *value, vector<double const *> const &args,
			 vector<vector<unsigned long> > const &dims) const
    {
	unsigned long nrow = NROW(dims[0]);
	unsigned long ncol = NCOL(dims[0]);
	unsigned long length = nrow * ncol;
	for (unsigned long i = 0; i < length; ++i) {
	    value[i] = args[0][(i / ncol) + (i % ncol) * nrow];
	}
    }

    bool Transpose::hasGradient(unsigned long i) const
    {
	return i == 0;
    }
    
    void Transpose::gradient(double *grad, vector<double const *> const &args,
			     vector<vector<unsigned long> > const &dims,
			     unsigned long i) const
    {
	/* A = B^T
	   where dimensions are:
	   A: P x Q
	   B: Q x P
	*/
	
	unsigned long Q = NROW(dims[0]);
	unsigned long P = NCOL(dims[0]);

	for (unsigned long i = 0; i < P; ++i) {
	    for (unsigned long j = 0; j < Q; ++j) {
		// Coordinates [i,j,k.l] with l==i, k==j
		grad[i + P*(j + Q*(j + Q*i))] += 1;
	    }
	}
    }
    
    vector<unsigned long> 
    Transpose::dim(vector <vector<unsigned long> > const &dims,
		   vector <double const *> const &) const
    {
	vector<unsigned long> ans(2);
	ans[0] = NCOL(dims[0]);
	ans[1] = NROW(dims[0]);
	return ans;
    }

    bool 
    Transpose::checkParameterDim (vector <vector<unsigned long> > const &dims) 
	const
    {
	return isScalar(dims[0]) || isVector(dims[0]) || isMatrix(dims[0]);
    }

    bool Transpose::isAdditive(vector<bool> const &, vector<bool> const &) const
    {
	return true;
    }
    
    bool Transpose::isScale(vector<bool> const &, vector<bool> const &) const
    {
	return true;
    }

    bool Transpose::isDiscreteValued(std::vector<bool> const &mask) const
    {
	return mask[0];
    }
}}
