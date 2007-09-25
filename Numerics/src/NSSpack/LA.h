#ifndef LA_H
#define LA_H

#include <config.h>

#define FCAST(T,X) (T *) (& X)
#define FCASTP(T,X) (T *) X

#ifdef HAVE_CBLAS_H
#include <cblas.h>
#ifndef HAVE_CLAPACK_H
#error "HAVE_CBLAS without HAVE_CLAPACK"
#endif
#include <clapack.h>

/* missing */
int clapack_dtrtrs(const enum ATLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum ATLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const enum CBLAS_DIAG Diag, const int n, const int nrhs, double *a, const int lda, double *b, const int ldb);

/* LA_ORDER, LA_SIDE specific cblas/clapack */
#ifndef LA_ORDER
#define LA_ORDER CblasColMajor
#endif
#ifndef LA_SIDE
#define LA_SIDE CblasLeft
#endif
/**/

#define LAPACK_NAME(N) clapack_##N
#define BLAS_NAME(N) cblas_##N
#define LA_TRANS CblasTrans
#define LA_NOTRANS CblasNoTrans
#define LA_COLMAJOR CblasColMajor
#define LA_ROWMAJOR CblasRowMajor
#define LA_UP CblasUpper
#define LA_LO CblasLower
#define LA_NONUNIT CblasNonUnit
#define LA_UNIT CblasUnit
#define INTEGER(N) N
#define INTEGERP(N) N
#define DOUBLE(X) X
#define T_UPLO(X) X
#define T_TRANS(X) X
#define T_DIAG(X) X
#define WITH_ORDER(ORDER, args...) \
  (ORDER, args)
#define LAPACK_4(F,A1,A2,A3,A4,INFO) \
  INFO = F(LA_ORDER,A1,A2,A3,A4)
#define LAPACK_7(F,A1,A2,A3,A4,A5,A6,A7,INFO) \
  INFO = F(LA_ORDER,A1,A2,A3,A4,A5,A6,A7)
#define LAPACK_9(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,INFO) \
  INFO = F(LA_ORDER,A1,A2,A3,A4,A5,A6,A7,A8,A9)
#define LAPACK_9_SIDED(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,INFO) \
  INFO = F(LA_ORDER,LA_SIDE,A1,A2,A3,A4,A5,A6,A7,A8,A9)


#else /* f2c or g2c + blaslapack.h */
#include "blaslapack.h"
#define BLAS_NAME(N) F77NAME(N)
#define LAPACK_NAME(N) F77NAME(N)
#define LA_TRANS "T"
#define LA_NOTRANS "N"
#define LA_COLMAJOR dummy
#define LA_ROWMAJOR dummy
#define LA_UP "U"
#define LA_LO "L"
#define LA_NONUNIT "N"
#define LA_UNIT "U"
#define INTEGER(N) FCAST(integer,N)
#define INTEGERP(N) FCASTP(integer,N)
#define DOUBLE(X) FCAST(double,X)
#define T_UPLO(X) FCASTP(char, X)
#define T_TRANS(X) FCASTP(char, X)
#define T_DIAG(X) FCASTP(char, X)
#define WITH_ORDER(ORDER, args...) (args)
#define APPLY(F,args...) F(args)
#define LAPACK_4 APPLY
#define LAPACK_7 APPLY
#define LAPACK_9 APPLY
#define LAPACK_9_SIDED APPLY
#endif /* HAVE_CBLAS */

/* DNRM2 - the euclidean norm of a vector via the function name, so
   that DNRM2 := sqrt( x'*x )
*/
#define DNRM2(N, X, INCX) \
  BLAS_NAME(dnrm2)(INTEGER(N), X, INTEGER(INCX))

/* DCOPY - a vector, x, to a vector, y
*/
#define DCOPY(N, X, INCX, Y, INCY) \
  BLAS_NAME(dcopy)(INTEGER(N), X, INTEGER(INCX), Y, INTEGER(INCY))

/* DGEMV - one of the matrix-vector operations y := alpha*A*x +
   beta*y, or y := alpha*A'*x + xbeta*y,
*/
#define DGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY) \
  BLAS_NAME(dgemv)WITH_ORDER(LA_ORDER, T_TRANS(TRANS), INTEGER(M), INTEGER(N), DOUBLE(ALPHA), A, INTEGER(LDA), X, INTEGER(INCX), DOUBLE(BETA), Y, INTEGER(INCY))

/*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',*/
#define DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) \
  BLAS_NAME(dgemm)WITH_ORDER(LA_ORDER, T_TRANS(TRANSA), T_TRANS(TRANSB), INTEGER(M), INTEGER(N), INTEGER(K), DOUBLE(ALPHA), A, INTEGER(LDA), B, INTEGER(LDB), DOUBLE(BETA), C, INTEGER(LDC))

/* DDOT - the dot product of two vectors
 */
#define DDOT(N, DX, INCX, DY, INCY) \
  BLAS_NAME(ddot)(INTEGER(N), DX, INTEGER(INCX), DY, INTEGER(INCY))

/* DAXPY - time a vector plus a vector
 */
#define DAXPY(N, DA, DX, INCX, DY, INCY) \
  BLAS_NAME(daxpy)(INTEGER(N), DOUBLE(DA), DX, INTEGER(INCX), DY, INTEGER(INCY))


/* DSCAL - a vector by a constant
 */
#define DSCAL(N,DA,DX,INCX) \
  BLAS_NAME(dscal)(INTEGER(N), DOUBLE(DA), DX, INTEGER(INCX))

/* DPOTRF - compute the Cholesky factorization of a real symmetric
   positive definite matrix A
 */
#define DPOTRF( UPLO, N, A, LDA, INFO ) \
  LAPACK_4(LAPACK_NAME(dpotrf), T_UPLO(UPLO), INTEGER(N), A , INTEGER(LDA), INTEGER(INFO))

/* DPOTRI - compute the inverse of a real symmetric positive definite
   matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
   computed by DPOTRF
*/
#define DPOTRI( UPLO, N, A, LDA, INFO) \
  LAPACK_4(LAPACK_NAME(dpotri), T_UPLO(UPLO), INTEGER(N), A, INTEGER(LDA), INTEGER(INFO))

/* DGESV - compute the solution to a real system of linear equations A * X = B,
 */
#define DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO ) \
  LAPACK_7(LAPACK_NAME(dgesv), INTEGER(N), INTEGER(NRHS), A, INTEGER(LDA), INTEGERP(IPIV), B, INTEGER(LDB), INTEGER(INFO))

/* DTRTRS - solve a triangular system of the form  A * X = B or A**T * X = B,
 */
#define DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO ) \
  LAPACK_9_SIDED(LAPACK_NAME(dtrtrs), T_UPLO(UPLO), T_TRANS(TRANS), T_DIAG(DIAG), INTEGER(N), INTEGER(NRHS), A, INTEGER(LDA), B, INTEGER(LDB), INTEGER(INFO) )

#endif /* LA_H */
