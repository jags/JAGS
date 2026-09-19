#include <config.h>

#include <string>
#include <cmath>
#include <algorithm>
#include <vector>

#include <module/ModuleError.h>
#include <util/integer.h>
#include <matrix/lapack.h>
#include <matrix/matrix.h>

using std::log;
using std::fabs;
using std::copy;
using std::vector;

namespace jags {

double logdet(double const *a, unsigned long n)
{
   // Log determinant of n x n symmetric positive-definite matrix a */
  
  unsigned long N = n*n;
  vector <double> acopy(N);
  copy(a, a + N, acopy.begin());

  vector<double> w(n);
  int lwork = -1;
  double worktest = 0;
  int info = 0;
  int ni = asInteger(n);
  jags_dsyev("N","L", &ni, &acopy[0], &ni, &w[0], &worktest, &lwork, &info);
  if (info != 0) {
    throwRuntimeError("unable to calculate workspace size for dsyev");
  }
  lwork = static_cast<int>(worktest);
  vector<double> work(static_cast<unsigned long>(lwork));
  jags_dsyev("N","L", &ni, &acopy[0], &ni, &w[0], &work[0], &lwork, &info);
  if (info != 0) {
    throwRuntimeError("unable to calculate eigenvalues in dsyev");
  }

  if (w[0] <= 0) {
      throwRuntimeError("Non positive definite matrix in call to logdet");
  }

  double logdet = 0;
  for (unsigned long i = 0; i < n; i++) {
    logdet += log(w[i]);
  }

  return logdet;
}

/*
bool check_symmetric_ispd(double const *a, unsigned long n)
{
    // Checks that an n x n symmetric matrix is positive definite.
    // The code is essentially the same as logdet, but we return
    //   false if the smallest eigenvalue is less than zero.
  
    unsigned long N = n*n;
    vector<double> acopy(N);
    copy(a, a+N, acopy.begin());

    //Workspace query to get optimal workspace
    vector<double> w(n);
    int lwork = -1;
    double worktest = 0;
    int info = 0;
    int ni = asInteger(n);
    F77_DSYEV("N","L", &ni, &acopy[0], &ni, &w[0], &worktest, &lwork, &info);
    if (info != 0) {
	throwRuntimeError("unable to calculate workspace size for dsyev");
    }
    lwork = static_cast<int>(worktest);
    vector<double> work(static_cast<unsigned long>(lwork));

    //Calculate eigenvalues
    F77_DSYEV("N","L", &ni, &acopy[0], &ni, &w[0], &work[0], &lwork, &info);
    if (info != 0) {
	throwRuntimeError("unable to calculate eigenvalues in dsyev");
    }

    return w[0] > 0;
}
*/

double det(double const *a, int n)
{
   // Determinant of n x n matrix a via the QR decomposition
  
  int N = n*n;
  vector<double> acopy(N);
  copy(a, a + N, acopy.begin());

  vector<double> tau(N);

  int lwork = -1;
  double worktest = 0;
  int info = 0;
  jags_dgeqrf(&n, &n, acopy.data(), &n, tau.data(), &worktest, &lwork, &info); 
  if (info != 0) {
    throwRuntimeError("unable to calculate workspace size for dgeqrf");
  }
  lwork = static_cast<int>(worktest);
  vector<double> work(lwork);
  jags_dgeqrf(&n, &n, acopy.data(), &n, tau.data(), work.data(), &lwork, &info); 
  if (info != 0) {
    throwRuntimeError("unable to calculate eigenvalues in dgeqrf");
  }

  double det = (n % 2 == 1) ? 1.0 : -1.0;
  for (int i = 0; i < n; i++) {
    det *= acopy[i + n*i];
  }

  return det;
}

  void cholesky (double *X, double const *A, unsigned long n)
  {
    unsigned long N = n*n;
    copy(A, A + N, X);
    
    int info = 0;
    int ni = asInteger(n);
    jags_dpotrf ("L", &ni, X, &ni, &info);
    if (info < 0) {
      throwLogicError("Illegal argument in cholesky");
    }
    else if (info > 0) {
      throwRuntimeError("Cannot calculate cholesky decomposition: not positive definite");
    }
  }

bool inverse_chol (double *X, double const *A, unsigned long n)
{
    /* invert n x n symmetric positive definite matrix A. Put result in X*/
    
    unsigned long N = n*n;
    copy(A, A + N, X);

    int info = 0;
    int ni = asInteger(n);
    jags_dpotrf ("L", &ni, X, &ni, &info);
    if (info != 0) return false;

    jags_dpotri ("L", &ni, X, &ni, &info); 
    if (info != 0) return false;

    //Copy lower to upper triangle
    for (unsigned long i = 0; i < n; ++i) {
	for (unsigned long j = 0; j < i; ++j) {
	    X[i*n + j] = X[j*n + i];
	}
    }

    return true;
}


bool inverse_lu (double *X, double const *A, unsigned long n)
{
    /* invert n x n matrix A. Put result in X */

    unsigned long N = n*n;
    vector<double> Acopy(N);
    copy(A, A + N, Acopy.begin());
    for (unsigned long i = 0; i < N; i++) {
	X[i] = 0;
    }
    for (unsigned long i = 0; i < n; i++) {
	X[i*n + i] = 1;
    }

    int info = 0;
    int ni = asInteger(n);
    vector<int> ipiv(n);
    jags_dgesv (&ni, &ni, Acopy.data(), &ni, ipiv.data(), X, &ni, &info);
    return info == 0;
}

    /*
bool check_symmetry(double const *x, unsigned long n, double tol)
{
    for (unsigned long i = 1; i < n; ++i) {
	double const *xp = x + i;
	double const *yp = x + n*i;
	for (unsigned long j = 0; j < i; ++j) {
	    if (fabs(*xp - *yp) > tol) return false;
	    xp += n;
	    yp++;
	}
    }
    return true;
}
    */

}
