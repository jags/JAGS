#include <config.h>
#include <sarray/util.h>
#include <sarray/SArray.h>
#include "matexp.h"
#include "DMState.h"

#include <cfloat>
#include <algorithm>
#include <cmath>

#include <Rmath.h>
#include <rng/RNG.h>

using std::min;
using std::max;
using std::vector;

#define INITIAL (*parameters[0]->value())
#define TIME (*parameters[1]->value())
#define NSTATE (parameters[2]->dim(true)[0])
#define INTENSITY (parameters[2]->value())

DMState::DMState()
  : DistDiscrete("dmstate", 3, DIST_SPECIAL, false) 
{}

DMState::~DMState()
{}

static void getParameters(int *initial, double *time,
			  int *nstates, double const **intensity, 
			  vector<SArray const *> const &parameters)
{
  *initial = static_cast<int>(*(parameters[0]->value()) + 1.0E-6);
  *time = *(parameters[1]->value());
  *nstates = parameters[2]->dim(true)[0];
  *intensity = parameters[2]->value();
}

bool DMState::checkParameterDim(vector<Index> const &dims) const
{
  return isScalar(dims[0]) && isScalar(dims[1]) && isSquareMatrix(dims[2]);
}

bool DMState::checkParameterDiscrete(vector<bool> const &mask) const
{
  return mask[0]; //Initial state must be discrete
}

bool 
DMState::checkParameterValue(vector<SArray const *> const &parameters) const
{
  int initial = 0;
  double time = 0;
  int nstates = 0;
  double const *intensity;
  getParameters(&initial, &time, &nstates, &intensity, parameters);

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

double DMState::l(vector<SArray const *> const &parameters,
		  bool absolute) const
{
    return 1;
}

double DMState::u(vector<SArray const *> const &parameters,
		  bool absolute) const
{
    return NSTATE;
}

/*
double DMState::logLikelihood(long x, vector<SArray const *> const &parameters) const
{
  int initial = 0;
  double time = 0;
  int nstate = 0;
  double const *intensity = 0;
  getParameters(&initial, &time, &nstate, &intensity, parameters);
  
  if (x < 1 || x > nstate) {
    return -DBL_MAX;
  }
  else {
    // FIXME: There is a big problem here with hidden Markov models.
    // If the node is not observed, then we will need to evaluate the
    // log-likelihood at different values of x, but with the same
    // parameters.  This is hugely inefficient at the moment because 
    // we calculate - and THROW AWAY - the transition probability matrix
    // each time.  

    //  SOLUTION: Maybe a cache?
    double *P = new double[nstate*nstate];
    //MatrixExp(P, intensity, nstate, time);
    MatrixExpPade(P, intensity, nstate, time);
    double lik = P[(initial - 1) + nstate * (x - 1)];
    delete [] P;
    if (lik <= 0) {
      // Allow for some numerical imprecision that may create small
      // negative probabilities
      return -DBL_MAX;
    }
    else {
      return log(lik);
    }
  }
}


long DMState::randomSample(vector<SArray const *> const &parameters) const
{
  int initial = 0;
  double time = 0;
  int nstate = 0;
  double const *intensity;
  getParameters(&initial, &time, &nstate, &intensity, parameters);
  
  double *P = new double[nstate*nstate];
  //MatrixExp(P, intensity, nstate, time);
  MatrixExpPade(P, intensity, nstate, time);

  // Draw the sample 
  double q = runif(0.0, 1.0);

  // Categorize 
  double sump = 0.0;
  for (long j = 1; j < nstate; j++) {
    sump += P[(initial - 1) + nstate * (j - 1)];
    if (q <= sump) {
      delete [] P;
      return j;
    }
  }

  delete [] P;
  return nstate;
}
*/

double DMState::d(double z, std::vector<SArray const *> const &parameters, 
	 bool give_log) const
{
  int initial = 0;
  double time = 0;
  int nstate = 0;
  double const *intensity = 0;
  getParameters(&initial, &time, &nstate, &intensity, parameters);
  
  long x = static_cast<long>(z);
  if (x < 1 || x > nstate) {
    return give_log ? -DBL_MAX : 0;
  }
  else {
    /*
      FIXME: There is a big problem here with hidden Markov models.
      If the node is not observed, then we will need to evaluate the
      log-likelihood at different values of x, but with the same
      parameters.  This is hugely inefficient at the moment because 
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
      return give_log ? -DBL_MAX : 0;
    }
    else {
      return give_log ? log(lik) : lik;
    }
  }
}

double DMState::p(double z, std::vector<SArray const *> const &parameters, 
		  bool lower, bool give_log) const
{
  int initial = 0;
  double time = 0;
  int nstate = 0;
  double const *intensity;
  getParameters(&initial, &time, &nstate, &intensity, parameters);
  
  long x = static_cast<long>(z);
  double sump = 0;
  if (x >= nstate) {
    sump = 1;
  }
  else {
    double *P = new double[nstate*nstate];
    MatrixExpPade(P, intensity, nstate, time);
    for (long j = 1; j <= x; ++j) {
      sump += P[(initial - 1) + nstate * (j - 1)];
    }
    delete [] P;
  }

  if (!lower)
    sump = max(1 - sump, 0.0);
  if (give_log)
    sump = log(sump);
  
  return sump; 
}

double DMState::q(double p, std::vector<SArray const *> const &parameters, bool lower,
                bool log_p) const
{
  if (log_p)
    p = exp(p);
  if (!lower)
    p = 1 - p;
    
  int initial = 0;
  double time = 0;
  int nstate = 0;
  double const *intensity;
  getParameters(&initial, &time, &nstate, &intensity, parameters);
  
  double *P = new double[nstate*nstate];
  MatrixExpPade(P, intensity, nstate, time);

  /* Categorize */
  double sump = 0.0;
  for (long j = 1; j < nstate; j++) {
    sump += P[(initial - 1) + nstate * (j - 1)];
    if (p <= sump) {
      delete [] P;
      return j;
    }
  }

  delete [] P;
  return nstate;
}

double DMState::r(std::vector<SArray const *> const &parameters, RNG *rng) 
  const
{
    return q(rng->uniform(), parameters, true, false);
}
