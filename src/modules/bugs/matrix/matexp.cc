/* 
   MatrixExp
   Function for calculating exponential of real square matrix.
   Based on pijt.c by Chris Jackson Copyright (C) 2001

   Changes:
   * Matrices are assumed to be filled by columns, not rows
   * Use LAPACK routines instead of EISPACK
   * const correctness
   * Ported to C++ with corresponding minor syntactic changes
*/ 

#include <config.h>

#include <stdexcept>
#include <cmath>

#include "lapack.h"
#include "matexp.h"

using std::pow;
using std::exp;
using std::log;
using std::runtime_error;

static int c_1 = 1;

/* index to treat a vector as a matrix, assuming matrix is filled by columns */
inline int MI(int i, int j, int ncols) {return j*ncols + i;}

/*
static void MatTranspose(double *AT, double const *A, int n)
{
  // Transpose a square matrix
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      AT[MI(i,j,n)] = A[MI(j,i,n)];
}
*/

static void MultMat(double *AB, double const *A, double const *B, int n)
{
  /* Multiply two square matrices together */
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      AB[MI(i, j, n)] = 0;
      for (int k = 0; k < n; ++k) 
	AB[MI(i, j, n)] += A[MI(i, k, n)] * B[MI(k, j, n)];
    }
  }
}

/*
static void 
MultMatDiag(double *AB, double const *diagA, double const *B, int n)
{
// Pre-multiply a square matrix by a diagonal matrix (given by a vector)
  for (int i = 0; i < (n*n); ++i) {
    AB[i] = 0;
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      AB[MI(i, j, n)] += diagA[i] * B[MI(i, j, n)];
    }
  }
}
*/

static void FormIdentity(double *A, int n)
{
  /* Set A to be an n x n identiy matrix */   
  for (int i = 0; i < (n*n); ++i) {
    A[i] = 0;
  }
  for (int i = 0; i < n; ++i) {
    A[MI(i, i, n)] = 1;
  }
}

/*
static bool repeated_entries(double const *vec, int n)
{
// Test whether a vector has any non-unique entries
  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      if (vec[j] == vec[i]) {
	return true;
      }
    }
  }
  return false;
}
*/

static void solve (double *X, double const *A, double const *B, int n)
{
  /* Solve AX = B, where all matrices are square */

  int N = n*n;
  double *Acopy = new double[N];
  F77_DCOPY(&N, A, &c_1, Acopy, &c_1);
  F77_DCOPY(&N, B, &c_1, X, &c_1);
  /*
  for (int i = 0; i < N; i++) {
    Acopy[i] = A[i];
    X[i] = B[i];
  }
  */
  int *ipiv = new int[n];
  int info = 0;
  F77_DGESV (&n, &n, Acopy, &n, ipiv, X, &n, &info);
  if (info != 0) {
    throw runtime_error("Unable to solve linear equations");
  }
  
  delete [] Acopy;
  delete [] ipiv;
}

static void MatrixExpSeries(double *expmat, double const *A, int n, double t)
{
  /* Calculate a matrix exponential using a power series approximation */
  /* Adapted from mexp in Jim Lindsey's rmutil library */  
  
  const int order = 20;   /* number of terms in series */
  const int underflow_correct = 3;
    
  double *At = new double[n*n];
  double *Apower = new double[n*n];
  double *Temp = new double[n*n];

  for (int i = 0; i < (n*n); ++i) {
    At[i] = A[i] * t / pow(2.0, underflow_correct);
  }
  FormIdentity(expmat, n);
  FormIdentity(Apower, n);
  for (int i = 1; i <= order; i++) {
    MultMat(Temp, At, Apower, n);
    for (int j = 0; j < (n*n); ++j){
      Apower[j] = Temp[j] / i;
      expmat[j] += Apower[j];
    }
  }
  for (int i = 0; i < underflow_correct; ++i){
    MultMat(Temp, expmat, expmat, n);
    for (int j = 0; j < (n*n); ++j)
      expmat[j] = Temp[j];
  }
    
  delete [] At;
  delete [] Apower;
  delete [] Temp;
}

/* Debugging: suppress the LAPACK routine */

void MatrixExp(double *expmat, double const *mat, int n, double t) {
    MatrixExpSeries(expmat, mat, n, t);
}

static void
padeseries (double *Sum, double const *A, int n, int order, 
            double scale, double *Temp)
{
  double zero = 0;
  int N = n*n;
  FormIdentity(Sum, n);
  for (int j = order; j >= 1; --j) {
    double s = (order-j+1) / (j*(2*order-j+1) * scale);
    F77_DGEMM ("n","n", &n, &n, &n, &s, Sum, &n, A, &n, &zero, Temp, &n);
    F77_DCOPY(&N, Temp, &c_1, Sum, &c_1);
    /*
    MultMat(Temp, Sum, A, n);
    for (int i = 0; i < N; ++i) {
      Sum[i] = Temp[i] * s;
    }
    */
    for (int r = 0; r < n; ++r) {
      Sum[r*n+r] += 1;
    }

  }
}

void 
MatrixExpPade(double *ExpAt, double const *A, const int n, double t)
{
  /* Calculate exp(A*t) by diagonal Pade approximation with scaling and
     squaring */

  int order = 8;
  int N = n*n;
  double *workspace = new double[4*N];
  double * const Temp = workspace;
  double * const At = workspace + N;
  double * const Num = workspace + 2*N;
  double * const Denom = workspace + 3*N;

  // Take a copy of the matrix and scale it by t
  F77_DCOPY(&N, A, &c_1, At, &c_1);
  F77_DSCAL(&N, &t, At, &c_1);
  /*
  for (int i = 0; i < N; ++i) {
    At[i] = A[i] * t;
  }
  */

  // Scale the matrix by a power of 2 
  double l1 = F77_DLANGE("1", &n, &n, At, &n, 0); //L-1 norm
  double linf = F77_DLANGE("i", &n, &n, At, &n, Temp); //L-Infinity norm
  /* 
     sqrt(l1 * linf) is an upper bound on the L2 norm of the matrix At
     (i.e the largest eigenvalue). We want to take the log, to base 2
     of this to get the smallest K, st ||At/2^K|| <= 1.
  */
  double K = (log(l1) + log(linf))/log(4.0);
  int npower = static_cast<int>(K) + 1;
  if (npower < 0) {
    npower = 0;
  }
  double scale = pow(2.0, npower);

  // Calculate exp(A/scale) by Pade series 

  padeseries (Num, At, n, order, scale, Temp);
  for (int i = 0; i < N; ++i) {
    At[i] = -At[i];
  }
  padeseries (Denom, At, n, order, scale, Temp);
  solve(ExpAt, Denom, Num, n);

  // Now repeatedly square the result

  for (int i = 0; i < npower; ++i) {
    /*
      F77_DCOPY(&N, ExpAt, &c_1, Temp, &c_1);
      double one = 1;
      double zero = 0;
      F77_DGEMM ("n","n", &n, &n, &n, &one, Temp, &n, Temp, &n, &zero, ExpAt,&n);
    */
    for (int j = 0; j < N; ++j) {
      Temp[j] = ExpAt[j];
    }
    MultMat(ExpAt, Temp, Temp, n);
  }

  delete [] workspace;
}
