#include <config.h>
#include <util/dim.h>

#include "Transpose.h"

using std::vector;

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
	unsigned long nrow = dims[0][0];
	unsigned long ncol = dims[0].size() == 2 ? dims[0][1] : 1;
	unsigned long length = nrow * ncol;
	for (unsigned long i = 0; i < length; ++i) {
	    value[i] = args[0][(i / ncol) + (i % ncol) * nrow];
	}
    }

    vector<unsigned long> 
    Transpose::dim (vector <vector<unsigned long> > const &dims,
		    vector <double const *> const &) const
    {
	vector<unsigned long> ans(2);
	ans[0] = dims[0].size() == 2 ? dims[0][1] : 1;
	ans[1] = dims[0][0];
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
