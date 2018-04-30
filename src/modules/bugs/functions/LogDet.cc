#include <config.h>

#include "LogDet.h"
#include "matrix.h"

#include <util/dim.h>
#include <util/integer.h>
#include <cmath>

using std::vector;

namespace jags {
namespace bugs {

    LogDet::LogDet ()
	: ArrayFunction ("logdet", 1)
    {
    }

    void LogDet::evaluate (double *x, vector<double const *> const &args,
			   vector<vector<unsigned long> > const &dims) const
    {
	*x = logdet(args[0], asInteger(dims[0][0]));
    }

    bool 
    LogDet::checkParameterDim (vector<vector<unsigned long> > const &dims) const
    {
	return isSquareMatrix(dims[0]) || isScalar(dims[0]);
    }

    vector<unsigned long>
    LogDet::dim(vector<vector<unsigned long> > const &dims,
		vector<double const *> const &values) const
    {
	return vector<unsigned long>(1,1);
    }

}}
