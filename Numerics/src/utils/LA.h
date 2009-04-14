#ifndef LA_H
#define LA_H

#include "NumericsConfig.h"

#define FCAST(T,X) (T *) (& X)
#define FCASTP(T,X) (T *) X

#ifndef FRAMEWORK_BLAS
#define OUTSIDE_FRAMEWORK_BLAS
#endif

#if defined(OUTSIDE_FRAMEWORK_BLAS) && defined(HAVE_CBLAS_H)
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
#define LAPACK_5(F,A1,A2,A3,A4,A5,INFO)   \
  INFO = F(LA_ORDER,A1,A2,A3,A4,A5)
#define LAPACK_6(F,A1,A2,A3,A4,A5,A6,INFO)  \
  INFO = F(LA_ORDER,A1,A2,A3,A4,A5,A6)
#define LAPACK_7(F,A1,A2,A3,A4,A5,A6,A7,INFO) \
  INFO = F(LA_ORDER,A1,A2,A3,A4,A5,A6,A7)
#define LAPACK_8(F,A1,A2,A3,A4,A5,A6,A7,A8,INFO)  \
  INFO = F(LA_ORDER,A1,A2,A3,A4,A5,A6,A7,A8)
#define LAPACK_9(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,INFO) \
  INFO = F(LA_ORDER,A1,A2,A3,A4,A5,A6,A7,A8,A9)
#define LAPACK_9_SIDED(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,INFO) \
  INFO = F(LA_ORDER,LA_SIDE,A1,A2,A3,A4,A5,A6,A7,A8,A9)
#define LAPACK_4_W LAPACK_4
#define LAPACK_5_W LAPACK_5

#else /* f2c or g2c + blaslapack.h */
#ifndef FRAMEWORK_BLAS
#include "blaslapack.h"
#endif
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
#define LAPACK_5 APPLY
#define LAPACK_6 APPLY
#define LAPACK_7 APPLY
#define LAPACK_8 APPLY
#define LAPACK_9 APPLY
#define LAPACK_9_SIDED APPLY
#define LAPACK_4_W(F,A1,A2,A3,A4,INFO) \
 ({ int C_LWORK=-1; \
    double * C_WORK; \
    C_WORK = malloc(sizeof *C_WORK); \
    if (C_WORK == NULL) -1; else { \
      F(A1,A2,A3,A4,C_WORK,INTEGER(C_LWORK),INFO); \
      C_LWORK = (int) (C_WORK[0]); \
      C_WORK = realloc(C_WORK, C_LWORK * sizeof *C_WORK); \
      F(A1,A2,A3,A4,C_WORK,INTEGER(C_LWORK),INFO); \
    } \
 })
#define LAPACK_5_W(F,A1,A2,A3,A4,A5,INFO) \
 ({ int C_LWORK=-1; \
    double * C_WORK; \
    C_WORK = malloc(sizeof *C_WORK); \
    if (C_WORK == NULL) -1; else { \
      F(A1,A2,A3,A4,A5,C_WORK,INTEGER(C_LWORK),INFO); \
      C_LWORK = (int) (C_WORK[0]); \
      C_WORK = realloc(C_WORK, C_LWORK * sizeof *C_WORK); \
      F(A1,A2,A3,A4,A5,C_WORK,INTEGER(C_LWORK),INFO); \
    } \
 })

#endif /* HAVE_CBLAS */

/* DNRM2 - the euclidean norm of a vector via the function name, so
   that DNRM2 := sqrt( x'*x )
*/
#define DNRM2(N, X, INCX) \
  ({ int C_N=N; \
     int C_INCX = INCX; \
     BLAS_NAME(dnrm2)(INTEGER(C_N), X, INTEGER(C_INCX)); \
  })

/* DCOPY - a vector, x, to a vector, y
*/
#define DCOPY(N, X, INCX, Y, INCY) \
  ({ int C_N=N; \
     int C_INCX=INCX; \
     int C_INCY=INCY; \
     BLAS_NAME(dcopy)(INTEGER(C_N), X, INTEGER(C_INCX), Y, INTEGER(C_INCY)); \
  })

/** DGEMV - one of the matrix-vector operations y := alpha*A*x +
   beta*y, or y := alpha*A'*x + xbeta*y,
   \param TRANS
*/
#define DGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY) \
  ({ int C_M = M; \
     int C_N = N; \
     double C_ALPHA = ALPHA; \
     int C_LDA = LDA; \
     int C_INCX = INCX; \
     double C_BETA = BETA; \
     int C_INCY = INCY; \
     BLAS_NAME(dgemv)WITH_ORDER(LA_ORDER, T_TRANS(TRANS), INTEGER(C_M), INTEGER(C_N), DOUBLE(C_ALPHA), A, INTEGER(C_LDA), X, INTEGER(C_INCX), DOUBLE(C_BETA), Y, INTEGER(C_INCY)); \
  })

/*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',*/
#define DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) \
  ({ int C_M = M; \
     int C_N = N; \
     int C_K = K; \
     double C_ALPHA = ALPHA; \
     int C_LDA = LDA; \
     int C_LDB = LDB; \
     double C_BETA = BETA; \
     int C_LDC = LDC; \
     BLAS_NAME(dgemm)WITH_ORDER(LA_ORDER, T_TRANS(TRANSA), T_TRANS(TRANSB), INTEGER(C_M), INTEGER(C_N), INTEGER(C_K), DOUBLE(C_ALPHA), A, INTEGER(C_LDA), B, INTEGER(C_LDB), DOUBLE(C_BETA), C, INTEGER(C_LDC)); \
  })

/* DDOT - the dot product of two vectors
 */
#define DDOT(N, DX, INCX, DY, INCY) \
  ({ int C_N = N; \
     int C_INCX = INCX; \
     int C_INCY = INCY; \
     BLAS_NAME(ddot)(INTEGER(C_N), DX, INTEGER(C_INCX), DY, INTEGER(C_INCY)); \
  })
/* DAXPY - time a vector plus a vector
 */
#define DAXPY(N, DA, DX, INCX, DY, INCY) \
  ({ int C_N = N; \
     double C_DA = DA; \
     int C_INCX = INCX; \
     int C_INCY = INCY; \
     BLAS_NAME(daxpy)(INTEGER(C_N), DOUBLE(C_DA), DX, INTEGER(C_INCX), DY, INTEGER(C_INCY)); \
  })


/* DSCAL - a vector by a constant
 */
#define DSCAL(N,DA,DX,INCX) \
  ({ int C_N = N; \
     double C_DA = DA; \
     int C_INCX = INCX; \
     BLAS_NAME(dscal)(INTEGER(C_N), DOUBLE(C_DA), DX, INTEGER(C_INCX)); \
  })

/* DPOTRF - compute the Cholesky factorization of a real symmetric
   positive definite matrix A
 */
#define DPOTRF( UPLO, N, A, LDA, INFO ) \
  ({ int C_N = N; \
     int C_LDA = LDA; \
     LAPACK_4(LAPACK_NAME(dpotrf), T_UPLO(UPLO), INTEGER(C_N), A , INTEGER(C_LDA), INTEGER(INFO)); \
  })

/* DPOTRI - compute the inverse of a real symmetric positive definite
   matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
   computed by DPOTRF
*/
#define DPOTRI( UPLO, N, A, LDA, INFO) \
  ({ int C_N = N; \
     int C_LDA = LDA; \
     LAPACK_4(LAPACK_NAME(dpotri), T_UPLO(UPLO), INTEGER(C_N), A, INTEGER(C_LDA), INTEGER(INFO)); \
  })

/* DGESV - compute the solution to a real system of linear equations A * X = B,
 */
#define DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO ) \
  ({ int C_N = N; \
     int C_NRHS = NRHS; \
     int C_LDA = LDA; \
     int C_LDB = LDB; \
     LAPACK_7(LAPACK_NAME(dgesv), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA), INTEGERP(IPIV), B, INTEGER(C_LDB), INTEGER(INFO)); \
  })

/* DGELS -
 *  DGELS solves overdetermined or underdetermined real linear systems
 *  involving an M-by-N matrix A, or its transpose, using a QR or LQ
 *  factorization of A.  It is assumed that A has full rank.
 */
#include "blaslapack.h"

#define DGELS( M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO  ) \
   ({int C_M = M; \
     int C_N = N; \
     int C_LWORK = LWORK; \
     double * C_WORK = WORK; \
     int C_NRHS = NRHS; \
     int C_LDA = LDA; \
     int C_LDB = LDB; \
     F77NAME(dgels)("N", FCAST(integer,C_M) , FCAST(integer,C_N), FCAST(integer,C_NRHS), FCASTP(double, A), FCAST(integer,C_LDA),FCASTP(double,B), FCAST(integer, C_LDB), FCASTP(double,C_WORK) ,FCAST(integer,LWORK), FCAST(integer,INFO)); \
  })

/* DTRTRS - solve a triangular system of the form  A * X = B or A**T * X = B,
 */
#define DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO ) \
  ({ int C_N = N; \
     int C_NRHS = NRHS; \
     int C_LDA = LDA; \
     int C_LDB = LDB; \
     LAPACK_9_SIDED(LAPACK_NAME(dtrtrs), T_UPLO(UPLO), T_TRANS(TRANS), T_DIAG(DIAG), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA), B, INTEGER(C_LDB), INTEGER(INFO) ); \
  })

/* DGETRF - LU factorization
 */
#define DGETRF(M,N,A,LDA,IPIV,INFO) \
  ({ int C_M = M; \
     int C_N = N; \
     int C_LDA = LDA; \
     LAPACK_5(LAPACK_NAME(dgetrf), INTEGER(C_M), INTEGER(C_N), A, INTEGER(C_LDA), INTEGERP(IPIV), INTEGER(INFO)); \
  })

/* DGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
 */
#define DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)  \
    ({ int C_N = N; \
     int C_NRHS = NRHS; \
     int C_LDA = LDA; \
     int C_LDB = LDB; \
     LAPACK_8(LAPACK_NAME(dgetrs), T_TRANS(TRANS), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA), INTEGERP(IPIV),B, INTEGER(C_LDB), INTEGER(INFO) ); \
  })

/* DGETRI - matrix inversion
 */
#define DGETRI(N,A,LDA,IPIV,INFO) \
  ({ int C_N = N; \
     int C_LDA = LDA; \
     LAPACK_4_W(LAPACK_NAME(dgetri), INTEGER(C_N), A, INTEGER(C_LDA), INTEGERP(IPIV),INTEGER(INFO)); \
  })


#endif /* LA_H */
