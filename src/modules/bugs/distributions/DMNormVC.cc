#include <config.h>
#include <util/dim.h>
#include <util/nainf.h>
#include "DMNormVC.h"
#include "DMNorm.h"

#include <lapack.h>
#include <matrix.h>

#include <cmath>
#include <vector>
#include <cfloat>

#include <JRmath.h>

#include "matrix.h"

using std::vector;

namespace jags {
    namespace bugs {
    
	DMNormVC::DMNormVC()
	    : ArrayDist("dmnorm.vcov", 2) 
	{}

	double
	DMNormVC::logDensity(double const *x, unsigned int m, PDFType type,
			     vector<double const *> const &parameters,
			     vector<vector<unsigned int> > const &dims,
			     double const *lower, double const *upper) const
	{
	    double const * mu = parameters[0];
	    double const * V  = parameters[1];

	    vector<double> T(m * m);
	    inverse_spd (&T[0], V, m);

	    double loglik = 0;
	    vector<double> delta(m);
	    for (unsigned int i = 0; i < m; ++i) {
		delta[i] = x[i] - mu[i];
		loglik -= (delta[i] * T[i + i * m] * delta[i])/2;
		for (unsigned int j = 0; j < i; ++j) {
		    loglik -= delta[i] * T[i + j * m] * delta[j];
		}
	    }

	    switch(type) {
	    case PDF_PRIOR:
		break;
	    case PDF_LIKELIHOOD:
		loglik -= logdet(V, m)/2;
		break;
	    case PDF_FULL:
		loglik -= logdet(V, m)/2 + m * M_LN_SQRT_2PI;
		break;
	    }
    
	    return loglik;
	}

	void
	DMNormVC::randomSample(double *x, unsigned int m,
			       vector<double const *> const &parameters,
			       vector<vector<unsigned int> > const &dims,
			       double const *lower, double const *upper,
			       RNG *rng) const
	{
	    double const * mu = parameters[0];
	    double const * T = parameters[1];
    
	    DMNorm::randomsample(x, mu, T, false, m, rng);
	}

	bool
	DMNormVC::checkParameterDim(vector<vector<unsigned int> > const &dims)
	    const
	{
	    //Allow scalar mean and precision. 
	    if (isScalar(dims[0]) && isScalar(dims[1])) return true;

	    //Vector mean and matrix precision
	    if (!isVector(dims[0])) return false;
	    if (!isSquareMatrix(dims[1])) return false;
	    if (dims[0][0] != dims[1][0]) return false;
    
	    return true;
	}

	vector<unsigned int>
	DMNormVC::dim(vector<vector<unsigned int> > const &dims) const
	{
	    return dims[0];
	}
	
	bool
	DMNormVC::checkParameterValue(vector<double const *> const &parameters,
				      vector<vector<unsigned int> > const &dims)
	    const
	{
	    double const *precision = parameters[1];
	    unsigned int n = dims[0][0];

	    return check_symmetry(precision, n) &&
		check_symmetric_ispd(precision, n);
	}


	void
	DMNormVC::support(double *lower, double *upper, unsigned int length,
			  vector<double const *> const &parameters,
			  vector<vector<unsigned int> > const &dims) const
	{
	    for (unsigned int i = 0; i < length; ++i) {
		lower[i] = JAGS_NEGINF;
		upper[i] = JAGS_POSINF;
	    }
	}

	void
	DMNormVC::typicalValue(double *x, unsigned int m,
			       vector<double const *> const &parameters,
			       vector<vector<unsigned int> > const &dims,
			       double const *lower, double const *upper) const
	{
	    for (unsigned int i = 0; i < m; ++i) {
		x[i] = parameters[0][i];
	    }
	}
	
	bool DMNormVC::isSupportFixed(vector<bool> const &fixmask) const
	{
	    return true;
	}

    }
}
