#include <config.h>

#include <stdexcept>
#include <cmath>

#include "lapack.h"
#include "matrix.h"

using std::log;
using std::runtime_error;

double logdet(double const *a, int n)
{
   // Log determinant of n x n symmetric positive matrix a */
  
  int N = n*n;
  double *acopy = new double[N];
  for (int i = 0; i < N; i++) {
    acopy[i] = a[i];
  }

  double *w = new double[n];
  int lwork = -1;
  double worktest = 0;
  int info = 0;
  F77_DSYEV("N","U", &n, acopy, &n, w, &worktest, &lwork, &info);
  if (info != 0) {
    delete [] acopy;
    delete [] w;
    throw runtime_error("unable to calculate workspace size for dsyev");
  }
  lwork = static_cast<int>(worktest);
  double *work = new double[lwork];
  F77_DSYEV("N","U", &n, acopy, &n, w, work, &lwork, &info);
  delete [] acopy;
  delete [] work;
  if (info != 0) {
    delete [] w;
    throw runtime_error("unable to calculate eigenvalues in dsyev");
  }

  if (w[0] <= 0) {
    throw runtime_error("Non positive definite matrix in call to logdet");
  }

  double logdet = 0;
  for (int i = 0; i < n; i++) {
    logdet += log(w[i]);
  }
  delete [] w;

  return logdet;
}

/*
double det(double const *a, int n)
{
   // Determinant of n x n matrix a via the QR decomposition
  
  int N = n*n;
  double *acopy = new double[N];
  for (int i = 0; i < N; i++) {
    acopy[i] = a[i];
  }

  int *jpvt = new int[n];
  for (int i = 0; i < n; i++) {
    jpvt[i] = 0;
  }

  double *tau = new double[N];

  int lwork = -1;
  double worktest = 0;
  int info = 0;
  F77_DGEQRF(&n, &n, acopy, &n, jpvt, tau, &worktest, &lwork, &info); 
  if (info != 0) {
    delete [] acopy;
    delete [] jpvt;
    delete [] tau;
    throw runtime_error("unable to calculate workspace size for dgeqrf");
  }
  lwork = static_cast<int>(worktest);
  double *work = new double[lwork];
  F77_DGEQRF(&n, &n, acopy, &n, jpvt, tau, work, &lwork, &info); 
  if (info != 0) {
    delete [] acopy;
    delete [] jpvt;
    delete [] tau;
    throw runtime_error("unable to calculate eigenvalues in dgeqrf");
  }

  double det = (n % 2 == 1) ? 1.0 : -1.0;
  for (int i = 0; i < n; i++) {
    det *= acopy[i + n*i];
  }

  delete [] acopy;
  delete [] jpvt;
  delete [] tau;

  return det;
}
*/


void inverse (double *X, double const *A, int n, bool spd)
{
  /* invert n x n smatrix A. Put result in X*/

  int N = n*n;
  double *Acopy = new double[N];
  for (int i = 0; i < N; i++) {
    Acopy[i] = A[i];
    X[i] = 0;
  }
  for (int i = 0; i < n; i++) {
    X[i*n + i] = 1;
  }

  int info = 0;
  if (spd) {
    //Symmetric, positive definite
    F77_DPOSV ("L", &n, &n, Acopy, &n, X, &n, &info); 
  }
  else {
    int *ipiv = new int[n];
    F77_DGESV (&n, &n, Acopy, &n, ipiv, X, &n, &info);
    delete [] ipiv;    
  }
  delete [] Acopy;
  if (info != 0) {
    throw runtime_error("Unable to invert matrix");
  }
}
