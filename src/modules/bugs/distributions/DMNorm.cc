#include <config.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <util/integer.h>

#include "DMNorm.h"

#include <lapack.h>
#include <matrix.h>

#include <cmath>
#include <vector>
#include <cfloat>

#include <JRmath.h>

using std::vector;

namespace jags {
namespace bugs {

DMNorm::DMNorm()
  : ArrayDist("dmnorm", 2) 
{}

double DMNorm::logDensity(double const *x, PDFType type,
			  vector<double const *> const &parameters,
			  vector<vector<unsigned long> > const &dims) const
{
    double const * mu = parameters[0];
    double const * T = parameters[1];
    unsigned long m = dims[0][0];
    
    double loglik = 0;
    vector<double> delta(m);
    for (unsigned long i = 0; i < m; ++i) {
	delta[i] = x[i] - mu[i];
	loglik -= (delta[i] * T[i + i * m] * delta[i])/2;
	for (unsigned long j = 0; j < i; ++j) {
	    loglik -= (delta[i] * T[i + j * m] * delta[j]);
	}
    }

    switch(type) {
    case PDF_PRIOR:
	break;
    case PDF_LIKELIHOOD:
	loglik += logdet(T, m)/2;
	break;
    case PDF_FULL:
	loglik += logdet(T, m)/2 - m * M_LN_SQRT_2PI;
	break;
    }
    
    return loglik;
}

void DMNorm::randomSample(double *x,
			  vector<double const *> const &parameters,
			  vector<vector<unsigned long> > const &dims,
			  RNG *rng) const
{
    double const * mu = parameters[0];
    double const * T = parameters[1];
    unsigned long m = dims[0][0];
    
    randomsample(x, mu, T, true, m, rng);
}

void DMNorm::randomsample(double *x, double const *mu, double const *T,
			  bool prec, unsigned long nrow, RNG *rng)
{
  //FIXME: do something with rng

  unsigned long N = nrow*nrow;
  vector<double> Tcopy(N);
  copy(T, T + N, Tcopy.begin());
  vector<double> w(nrow);

  int info = 0;
  double worktest;
  int lwork = -1;
  int nr = asInteger(nrow);
  // Workspace query
  F77_DSYEV ("V", "L", &nr, &Tcopy[0], &nr, &w[0], &worktest, &lwork, &info);
  // Now get eigenvalues/vectors with optimal work space
  lwork = static_cast<int>(worktest);
  vector<double> work(static_cast<unsigned long>(lwork));
  F77_DSYEV ("V", "L", &nr, &Tcopy[0], &nr, &w[0], &work[0], &lwork, &info);

  /* Generate independent random normal variates, scaled by
     the eigen values. We reuse the array w. */
  if (prec) {
      for (unsigned long i = 0; i < nrow; ++i) {
	  w[i] = rnorm(0, 1/sqrt(w[i]), rng);
      }
  }
  else {
      for (unsigned long i = 0; i < nrow; ++i) {
	  w[i] = rnorm(0, sqrt(w[i]), rng);
      }
  }

  /* Now transform them to dependant variates 
    (On exit from DSYEV, Tcopy contains the eigenvectors)
  */
  for (unsigned long i = 0; i < nrow; ++i) {
      x[i] = mu ? mu[i] : 0;
      for (unsigned long j = 0; j < nrow; ++j) {
	  x[i] += Tcopy[i + j * nrow] * w[j]; 
      }
  }
}

bool DMNorm::checkParameterDim(vector<vector<unsigned long> > const &dims) const
{
    //Allow scalar mean and precision. 
    if (isScalar(dims[0]) && isScalar(dims[1])) return true;

    //Vector mean and matrix precision
    if (!isVector(dims[0])) return false;
    if (!isSquareMatrix(dims[1])) return false;
    if (dims[0][0] != dims[1][0]) return false;
    
    return true;
}

bool DMNorm::checkParameterValue(vector<double const *> const &,
				 vector<vector<unsigned long> > const &)
    const
{
    return true; //FIXME: Define default in base clase
}
    
vector<unsigned long> DMNorm::dim(vector<vector<unsigned long> > const &dims) const
{
    return dims[0];
}

void DMNorm::support(double *lower, double *upper,
		     vector<double const *> const &,
		     vector<vector<unsigned long> > const &dims) const
{
    unsigned long length = dims[0][0];
    for (unsigned long i = 0; i < length; ++i) {
	lower[i] = JAGS_NEGINF;
	upper[i] = JAGS_POSINF;
    }
}

bool DMNorm::isSupportFixed(vector<bool> const &) const
{
    return true;
}

}}
