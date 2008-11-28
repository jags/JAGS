#include <config.h>
#include <util/dim.h>
#include <util/nainf.h>
#include "DMNorm.h"

#include <lapack.h>
#include <matrix.h>

#include <cmath>
#include <vector>
#include <stdexcept>

#include <JRmath.h>

using std::vector;
using std::logic_error;

DMNorm::DMNorm()
  : Distribution("dmnorm", 2, false, false) 
{}

double DMNorm::logLikelihood(double const *x, unsigned int m,
			     vector<double const *> const &parameters,
			     vector<vector<unsigned int> > const &dims,
			     double const *lower, double const *upper) const
{
    double const * mu = parameters[0];
    double const * T = parameters[1];

    double loglik = logdet(T, m)/2;
    double * delta = new double[m];
    for (unsigned int i = 0; i < m; ++i) {
	delta[i] = x[i] - mu[i];
	loglik -= (delta[i] * T[i + i * m] * delta[i])/2;
	for (unsigned int j = 0; j < i; ++j) {
	    loglik -= (delta[i] * T[i + j * m] * delta[j]);
	}
    }
    delete [] delta;

    return loglik;
}

void DMNorm::randomSample(double *x, unsigned int length,
			  vector<double const *> const &parameters,
			  vector<vector<unsigned int> > const &dims,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
    double const * mu = parameters[0];
    double const * T = parameters[1];
    
    randomsample(x, mu, T, true, length, rng);
}

void DMNorm::randomsample(double *x, double const *mu, double const *T,
			  bool prec, int nrow, RNG *rng)
{
  //FIXME: do something with rng

  int N = nrow*nrow;
  double * Tcopy = new double[N];
  for (int i = 0; i < N; ++i) {
    Tcopy[i] = T[i];
  }
  double * w = new double[nrow];

  int info = 0;
  double worktest;
  int lwork = -1;
  // Workspace query
  F77_DSYEV ("V", "L", &nrow, Tcopy, &nrow, w, &worktest, &lwork, &info);
  // Now get eigenvalues/vectors with optimal work space
  lwork = static_cast<int>(worktest + DBL_EPSILON);
  double * work = new double[lwork];
  F77_DSYEV ("V", "L", &nrow, Tcopy, &nrow, w, work, &lwork, &info);
  delete [] work;

  /* Generate independent random normal variates, scaled by
     the eigen values. We reuse the array w. */
  if (prec) {
      for (int i = 0; i < nrow; ++i) {
	  w[i] = rnorm(0, 1/sqrt(w[i]), rng);
      }
  }
  else {
      for (int i = 0; i < nrow; ++i) {
	  w[i] = rnorm(0, sqrt(w[i]), rng);
      }
  }

  /* Now transform them to dependant variates 
    (On exit from DSYEV, Tcopy contains the eigenvectors)
  */
  for (int i = 0; i < nrow; ++i) {
      x[i] = mu ? mu[i] : 0;
      for (int j = 0; j < nrow; ++j) {
	  x[i] += Tcopy[i + j * nrow] * w[j]; 
      }
  }
  delete [] w;
  delete [] Tcopy;
}

bool DMNorm::checkParameterDim(vector<vector<unsigned int> > const &dims) const
{
  if (!isVector(dims[0]))
    return false;
  if (!isSquareMatrix(dims[1]))
    return false;
  if (dims[0][0] != dims[1][0])
    return false;
 
  return true;
}

vector<unsigned int> DMNorm::dim(vector<vector<unsigned int> > const &dims) const
{
    return dims[0];
}

bool
DMNorm::checkParameterValue(vector<double const *> const &parameters,
			    vector<vector<unsigned int> > const &dims) const
{
    unsigned int n = dims[0][0];
    double const *T = parameters[1];
    // Check symmetry
    for (unsigned int i = 1; i < n; i++) {
	for (unsigned int j = 0; j < i - 1; j++) {
	    if (fabs(T[i + j*n] - T[j + i*n]) > 16 * DBL_EPSILON)
		return false;
	}
    }
    // Don't bother checking positive definiteness

    return true;
}


void DMNorm::support(double *lower, double *upper, unsigned int length,
		     vector<double const *> const &parameters,
		     vector<vector<unsigned int> > const &dims) const
{
    for (unsigned int i = 0; i < length; ++i) {
	lower[i] = JAGS_NEGINF;
	upper[i] = JAGS_POSINF;
    }
}

void DMNorm::typicalValue(double *x, unsigned int length,
			  vector<double const *> const &parameters,
			  vector<vector<unsigned int> > const &dims,
			  double const *lower, double const *upper) const
{
    for (unsigned int i = 0; i < length; ++i) {
	x[i] = parameters[0][i];
    }
}

bool DMNorm::isSupportFixed(vector<bool> const &fixmask) const
{
    return true;
}
