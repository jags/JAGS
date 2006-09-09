#ifndef JAGS_LAPACK_H_
#define JAGS_LAPACK_H_

#define F77_DGEEV F77_FUNC(dgeev,DGEEV)
#define F77_DSYEV F77_FUNC(dsyev,DSYEV)
#define F77_DGESV F77_FUNC(dgesv,DGESV)
#define F77_DSYSV F77_FUNC(dsysv,DSYSV)
#define F77_DPOSV F77_FUNC(dposv,DPOSV)
/* #define F77_DGEQP3 F77_FUNC(dgeqp3,DGEQP3) */
#define F77_DPOTRF F77_FUNC(dpotrf,DPOTRF)
#define F77_DLANGE F77_FUNC(dlange,DLANGE)
#define F77_DGEMM  F77_FUNC(dgemm,DGEMM)
#define F77_DCOPY  F77_FUNC(dcopy,DCOPY)
#define F77_DSCAL  F77_FUNC(dscal,DSCAL)

extern "C" {
	void F77_DGEEV (const char* jobvl, const char* jobvr,
			const int* n, double* a, const int* lda,
			double* wr, double* wi, double* vl, const int* ldvl,
			double* vr, const int* ldvr,
			double* work, const int* lwork, int* info);

	void F77_DSYEV (const char* jobz, const char* uplo,
			const int* n, double* a, const int* lda,
			double* w, 
			double* work, const int* lwork, int* info);
	
	
	void F77_DGESV (const int* n, const int* nrhs, double* a, 
			const int* lda,	int* ipiv, double* b, const int* ldb,
			int* info);

	void F77_DSYSV (const char *uplo, const int *n, const int *nrhs,
			double *a, const int *lda, int *ipiv, double *b,
			const int *ldb, double *work, const int *lwork,
			int *info);
	
	void F77_DPOSV (const char* uplo, const int *n, const int *nrhs,
			double* a, const int *lda, double* b, const int* ldb,
			int* info);

	/*
	void F77_DGEQP3 (const int* m, const int* n, double* a,
			 const int *lda, int *jpvt, double *tau,
			 double *work, const int *lwork, int *info);

	void F77_DGEQRF (const int *m, const *int n, double *a,
			 const int *lda, double *tau, double *work,
			 const int *lwork);
	*/

	void F77_DPOTRF (const char *uplo, const int *n, double *a,
			 const int *lda, const int *info);

  double F77_DLANGE (const char *norm, const int *m, const int *n,
		     const double *a, const int *lda, double *work);


  /* BLAS routines */

  void F77_DGEMM (const char* transa, const char* transb, const int* m,
		  const int* n , const int* k, const double* alpha,
		  const double* A, const int* lda, 
		  const double* B, const int* ldb,
		  const double* beta, double* C, const int* ldc);

  void F77_DCOPY (const int* n, double const *X, const int* incx, double *Y,
                  const int *incy);

  void F77_DSCAL (const int* n, double const *alpha, double *X, 
                  const int *incx);
}

#endif
