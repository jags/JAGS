#include <config.h>
#include <lapack.h>
#include <matrix.h>
#include <sarray/util.h>
#include <sarray/SArray.h>
#include "DMNorm.h"

#include <cmath>
#include <vector>
#include <stdexcept>

#include <Rmath.h>

using std::vector;
using std::logic_error;

DMNorm::DMNorm()
  : Distribution("dmnorm", 2, false, false) 
{}

DMNorm::~DMNorm()
{}

double DMNorm::logLikelihood(SArray const &x,
			     vector<SArray const *> const &parameters) const
{
  double const *y = x.value();
  int m = parameters[0]->length();
  double const * mu = parameters[0]->value();
  double const * T = parameters[1]->value();

  double loglik = logdet(T, m)/2;
  double * delta = new double[m];
  for (int i = 0; i < m; ++i) {
    delta[i] = y[i] - mu[i];
    loglik -= (delta[i] * T[i + i * m] * delta[i])/2;
    for (int j = 0; j < i; ++j) {
      loglik -= (delta[i] * T[i + j * m] * delta[j]);
    }
  }
  delete [] delta;

  return loglik;
}

void DMNorm::randomSample(SArray &x,
			  vector<SArray const *> const &parameters,
			  RNG *rng) const
{
    double const * mu = parameters[0]->value();
    double const * T = parameters[1]->value();
    int nrow = parameters[0]->length();

    double *y = new double[nrow];
    randomsample(y, mu, T, nrow, rng);
    x.setValue(y, nrow);
    delete [] y;
}

void DMNorm::randomsample(double *x, double const *mu, double const *T,
			  int nrow, RNG *rng)
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
  for (int i = 0; i < nrow; ++i) {
    w[i] = rnorm(0, 1/sqrt(w[i]), rng);
  }

  /* Now transform them to dependant variates 
    (On exit from DSYEV, Tcopy contains the eigenvectors)
  */
  for (int i = 0; i < nrow; ++i) {
    x[i] = mu[i];
    for (int j = 0; j < nrow; ++j) {
      x[i] += Tcopy[i + j * nrow] * w[j];
    }
  }
  delete [] w;
  delete [] Tcopy;
}

bool DMNorm::checkParameterDim(vector<Index> const &dims) const
{
  if (!isVector(dims[0]))
    return false;
  if (!isSquareMatrix(dims[1]))
    return false;
  if (dims[0][0] != dims[1][0])
    return false;
 
  return true;
}

Index DMNorm::dim(vector<Index> const &dims) const
{
  return dims[0];
}

bool DMNorm::checkParameterValue(vector<SArray const *> const &parameters) const
{
  long n = parameters[0]->length();

  double const *T = parameters[1]->value();
  // Check symmetry
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      if (fabs(T[i + j*n] - T[j + i*n]) > DBL_EPSILON)
	return false;
    }
  }
  // Don't bother checking positive definiteness

  return true;
}

unsigned long DMNorm::df(std::vector<SArray const *> const &parameters) const
{
    return parameters[0]->length();
}

double 
DMNorm::lowerSupport(unsigned long i,
		     std::vector<SArray const *> const &parameters) const
{
  unsigned long m = parameters[0]->length();
  if (i >= m)
    throw logic_error("Invalid index in DMNorm::lowerSupport");
  
  return -DBL_MAX;
}

double 
DMNorm::upperSupport(unsigned long i,
		     std::vector<SArray const *> const &parameters) const
{
  unsigned long m = parameters[0]->length();
  if (i >= m)
    throw logic_error("Invalid index in DMNorm::upperSupport");

  return DBL_MAX;
}

void 
DMNorm::typicalValue(SArray &x, std::vector<SArray const *> const &parameters)
  const
{
  x.setValue(parameters[0]->value(), parameters[0]->length());
}
