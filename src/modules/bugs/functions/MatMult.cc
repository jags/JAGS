#include <config.h>

#include "MatMult.h"
#include <util/dim.h>
#include <util/integer.h>

#include "lapack.h"

using std::vector;

namespace jags {
namespace bugs {

    //FIXME: deparse?

    MatMult::MatMult()
	: ArrayFunction("%*%", 2)
    {
    }

    void 
    MatMult::evaluate (double *value, vector<double const *> const &args,
		       vector<vector<unsigned long> > const &dims) const
    {
	int d1, d2, d3;

	if (dims[0].size() == 1) {
	    d1 = 1;
	    d2 = asInteger(dims[0][0]);
	}
	else {
	    d1 = asInteger(dims[0][0]);
	    d2 = asInteger(dims[0][1]);
	}
	if (dims[1].size() == 1) {
	    d3 = 1;
	}
	else {
	    d3 = asInteger(dims[1][1]);
	}
    
	double one = 1, zero = 0;
	F77_DGEMM ("N", "N", &d1, &d3, &d2, &one,
		   args[0], &d1, args[1], &d2, &zero, value, &d1);
    }

    vector<unsigned long> 
    MatMult::dim (vector <vector<unsigned long> > const &dims,
		  vector<double const *> const &) const
    {
	vector<unsigned long> ans(2,1);

	if (dims[0].size() == 2) {
	    ans[0] = dims[0][0];
	}
	if (dims[1].size() == 2) {
	    ans[1] = dims[1][1];
	}

	return drop(ans);
    }

    bool 
    MatMult::checkParameterDim (vector<vector<unsigned long> > const &dims) const
    {
	if (dims[0].size() > 2 || dims[1].size() > 2) {
	    return false;
	}
	
	if (dims[0].size() == 1) {
	    return dims[0][0] == dims[1][0];
	}
	else {
	    return dims[0][1] == dims[1][0];
	}
    }


    bool 
    MatMult::isScale(vector<bool> const &mask, vector<bool> const &fix) const
    {
	//Test for quadratic terms
	if (mask[0] && mask[1]) {
	    return false;
	}
    
	if (fix.empty()) {
	    return true;
	}
	else {
	    return (mask[0] || fix[0]) && (mask[1] || fix[1]);
	}
    }

    bool MatMult::isDiscreteValued(vector<bool> const &mask) const
    {
	return mask[0] && mask[1];
    }
}}
