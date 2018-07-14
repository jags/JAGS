#include <config.h>
#include "DMState.h"
#include <util/dim.h>
#include <util/nainf.h>
#include "matexp.h"


#include <cfloat>
#include <algorithm>
#include <cmath>

#include <JRmath.h>
#include <rng/RNG.h>

using std::min;
using std::max;
using std::vector;

#define INITIAL(par)   static_cast<int>(*par[0])
#define TIME(par)      (*par[1])
#define NSTATE(dims)   (dims[2][0])
#define INTENSITY(par) (par[2])

namespace jags {
namespace msm {

DMState::DMState()
    : ArrayDist ("dmstate", 3) 
{}

bool 
DMState::checkParameterDim(vector<vector<unsigned long> > const &dims) const
{
    return isScalar(dims[0]) && isScalar(dims[1]) && isSquareMatrix(dims[2]);
}

bool DMState::checkParameterDiscrete(vector<bool> const &mask) const
{
    return mask[0]; //Initial state must be discrete
}

bool 
DMState::checkParameterValue(vector<double const *> const &par,
			     vector<vector<unsigned long> > const &dims) const
{
    int initial = static_cast<int>(INITIAL(par));
    double time = TIME(par);
    int nstates = NSTATE(dims);
    double const *intensity = INTENSITY(par);

    // Number of states must be greater than 1
    if (nstates <= 1)
	return false;

    // Current state must be in range
    if (initial < 1 || initial > nstates)
	return false;

    // Time must be positive
    if (time < 0)
	return false;

    // Intensity parameter must be a real intensity matrix
    for (int i = 0; i < nstates; ++i) {
	double sum = 0.0;
	for (int j = 0; j < nstates; ++j) {
	    double rate = intensity[i + nstates * j];
	    // Diagonal elements are non-positive. Others are non-negative
	    if (i == j) {
		if (rate > 0.0)
		    return false;
	    }
	    else {
		if (rate < 0.0)
		    return false;
	    }
	    sum += rate;
	}
	// Rows of intensity matrix sum to 0
	if (fabs(sum) > 1.0E-6) {
	    // fixme: we should relate the tolerance to the machine precision
	    return false;
	}
    }

    return true;
}


double DMState::logDensity(double const *y, PDFType type,
		    vector<double const *> const &par,
		    vector<vector<unsigned long> > const &dims,
		    double const *lower, double const *upper) const
{
    int initial = INITIAL(par);
    double time = TIME(par);
    int nstate = NSTATE(dims);
    double const *intensity = INTENSITY(par);
  
    int x = static_cast<int>(*y);
    if (x < 1 || x > nstate) {
	return JAGS_NEGINF;
    }
    else {
	/*
	  FIXME: There is a big problem here with hidden Markov models.
	  If the node is not observed, then we will need to evaluate the
	  log-likelihood at different values of x, but with the same
	  par.  This is hugely inefficient at the moment because 
	  we calculate - and THROW AWAY - the transition probability matrix
	  each time.  

	  SOLUTION: Maybe a cache?
	*/
	double *P = new double[nstate*nstate];
	//MatrixExp(P, intensity, nstate, time);
	MatrixExpPade(P, intensity, nstate, time);
	double lik = P[(initial - 1) + nstate * (x - 1)];
	delete [] P;
	if (lik <= 0) {
	    /*
	      Allow for some numerical imprecision that may create small
	      negative probabilities
	    */
	    return JAGS_NEGINF;
	}
	else {
	    return log(lik);
	}
    }
}

static double q(double p, int initial, double time, unsigned int nstate,
		double const *intensity)
{
    double *P = new double[nstate*nstate];
    MatrixExpPade(P, intensity, nstate, time);
    
    /* Categorize */
    double sump = 0.0;
    for (unsigned int j = 1; j < nstate; j++) {
	sump += P[(initial - 1) + nstate * (j - 1)];
	if (p <= sump) {
	    delete [] P;
	    return j;
	}
    }

    delete [] P;
    return nstate;
}

void DMState::randomSample(double *x,
			   vector<double const *> const &par,
			   vector<vector<unsigned long> > const &dims,
			   double const *lower, double const *upper,
			   RNG *rng) const
{
    *x = q(runif(0, 1, rng), INITIAL(par), TIME(par), NSTATE(dims),
	   INTENSITY(par));
}

bool DMState::isSupportFixed(vector<bool> const &fixmask) const
{
    return true;
}

void DMState::support(double *lower, double *upper,
		      vector<double const *> const &par,
		      vector<vector<unsigned long> > const &dims) const
{
    *lower = 1;
    *upper = NSTATE(dims);
}

    vector<unsigned long> 
    DMState::dim(vector<vector<unsigned long> > const &dims) const
    {
	return vector<unsigned long>(1, 1);
    }


bool DMState::isDiscreteValued(vector<bool> const &mask) const
{
    return true;
}

}}
