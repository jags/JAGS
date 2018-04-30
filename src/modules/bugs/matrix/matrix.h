#ifndef MATRIX_H_
#define MATRIX_H_

namespace jags {
namespace bugs {

/**
 * Inverts a general square matrix using the LAPACK routine DGESV
 *
 * @param X Pointer to an array of length n squared, which will contain
 * the inverse on exit.
 *
 * @param A pointer to array containing the values of the matrix
 *
 * @param n number or rows or columns in the matrix
 */
bool inverse_lu (double *X, double const *A, unsigned long n);

/**
 * Inverts a symmetrix positive definite matrix by Cholesky
 * decomposition using the LAPACK routines DPOTRF and DPOTRI.
 * 
 * @param X Pointer to an array of length n squared, which will contain
 * the inverse on exit.
 *
 * @param A pointer to array containing the values of the matrix. Only
 * the lower triangle of the matrix (in column-major order) is used.
 *
 * @param n number or rows or columns in the matrix
 */
bool inverse_chol (double *X, double const *A, unsigned long n);

/**
 * Log determinant of a symmetric positive definite matrix
 *
 * @param A pointer to array containing the values of the matrix. Only
 * the lower triangle (in column-major order) is used.
 *
 * @param n number or rows or columns in the matrix
 */
double logdet(double const *A, unsigned long n);

}}

#endif /* MATRIX_H_ */
