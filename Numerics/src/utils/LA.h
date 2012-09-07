/* Siconos-Numerics, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#ifndef LA_H
#define LA_H

#include <assert.h>
#include "NumericsConfig.h"
#include <stdlib.h>

#if defined(USE_MKL)
#include<SiconosMKL.h>

#else /* do not use MKL */

#define FCAST(T,X) (T *) (& X)
#define FCASTP(T,X) (T *) X

#ifndef FRAMEWORK_BLAS
#define OUTSIDE_FRAMEWORK_BLAS
#endif

#if defined(OUTSIDE_FRAMEWORK_BLAS) && defined(HAVE_CBLAS_H) && defined(HAVE_CLAPACK_H)


#if defined(__cplusplus) && !defined (_NUMERICS_INTERNAL_CXX_)
#if defined(HAVE_ATLAS)
extern "C"
{
#include <atlas/cblas.h>
#include <atlas/clapack.h>
}
#else /* HAVE_ATLAS */
#include <cblas.h>
#include <clapack.h>
#endif /* HAVE_ATLAS */

#else /* __cplusplus */
#include <cblas.h>
#include <clapack.h>
#endif /* __cplusplus */

/* missing */
#ifdef HAVE_ATLAS
int clapack_dtrtrs(const enum ATLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum ATLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const enum CBLAS_DIAG Diag, const int n, const int nrhs, double *a, const int lda, double *b, const int ldb);
#endif

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
#define CLA_TRANS CblasTrans
#define CLA_NOTRANS CblasNoTrans
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
#define WITH_ORDER(ORDER, ...) \
  (ORDER, __VA_ARGS__)
#define LAPACK_4(F,A1,A2,A3,A4,INFO) \
  INFO = F(LA_ORDER,A1,A2,A3,A4)
#define LAPACK_5(F,A1,A2,A3,A4,A5,INFO)  \
  INFO = F(LA_ORDER,A1,A2,A3,A4,A5)
#define LAPACK_6(F,A1,A2,A3,A4,A5,A6,INFO) \
  INFO = F(LA_ORDER,A1,A2,A3,A4,A5,A6)
#define LAPACK_7(F,A1,A2,A3,A4,A5,A6,A7,INFO) \
  INFO = F(LA_ORDER,A1,A2,A3,A4,A5,A6,A7)
#define LAPACK_8(F,A1,A2,A3,A4,A5,A6,A7,A8,INFO) \
  INFO = F(LA_ORDER,A1,A2,A3,A4,A5,A6,A7,A8)
#define LAPACK_9(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,INFO) \
  INFO = F(LA_ORDER,A1,A2,A3,A4,A5,A6,A7,A8,A9)
#define LAPACK_9_SIDED(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,INFO) \
  INFO = F(LA_ORDER,LA_SIDE,A1,A2,A3,A4,A5,A6,A7,A8,A9)
#define LAPACK_4_W LAPACK_4
#define LAPACK_5_W LAPACK_5

#else /* f2c or g2c + blaslapack.h */
#include "blaslapack.h"
#define BLAS_NAME(N) F77NAME(N)
#define LAPACK_NAME(N) F77NAME(N)
#define LA_TRANS "T"
#define LA_NOTRANS "N"
#define CLA_TRANS 'T'
#define CLA_NOTRANS 'N'
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
#define WITH_ORDER(ORDER, ...) (__VA_ARGS__)
#define APPLY(F, ...) F(__VA_ARGS__)
#define LAPACK_4 APPLY
#define LAPACK_5 APPLY
#define LAPACK_6 APPLY
#define LAPACK_7 APPLY
#define LAPACK_8 APPLY
#define LAPACK_9 APPLY
#define LAPACK_9_SIDED APPLY
#define LAPACK_4_W(F,A1,A2,A3,A4,INFO)                    \
  int C_LWORK=-1;                                      \
double* C_WORK;                                      \
C_WORK = (double*)malloc(sizeof *C_WORK);                      \
assert(C_WORK);                                       \
F(A1,A2,A3,A4,C_WORK,INTEGER(C_LWORK),INFO);          \
C_LWORK = (int) (C_WORK[0]);                          \
C_WORK = (double*)realloc(C_WORK, C_LWORK * sizeof *C_WORK);   \
F(A1,A2,A3,A4,C_WORK,INTEGER(C_LWORK),INFO);          \
free(C_WORK);                                         \
 
#define LAPACK_5_W(F,A1,A2,A3,A4,A5,INFO) \
  ({ int C_LWORK=-1;                                      \
   double * C_WORK;                                      \
   C_WORK = malloc(sizeof *C_WORK);                      \
   assert(C_WORK);                                       \
   F(A1,A2,A3,A4,A5,C_WORK,INTEGER(C_LWORK),INFO);       \
   C_LWORK = (int) (C_WORK[0]);                          \
   C_WORK = realloc(C_WORK, C_LWORK * sizeof *C_WORK);   \
   F(A1,A2,A3,A4,A5,C_WORK,INTEGER(C_LWORK),INFO);       \
   free(C_WORK);                                         \
   })

#endif /* HAVE_CBLAS */

/* DNRM2 - the euclidean norm of a vector via the function name, so
   that DNRM2 := sqrt( x'*x )
   */
static inline double DNRM2(int N, double* X, int INCX)
{
  int C_N = N;
  int C_INCX = INCX;
  assert(C_N > 0);
  assert(C_INCX > 0);
  assert(X != NULL);
  return BLAS_NAME(dnrm2)(INTEGER(C_N), X, INTEGER(C_INCX));
}


/* DCOPY - a vector, x, to a vector, y
*/
static inline void DCOPY(int N, double* X, int INCX, double* Y, int INCY)
{
  int C_N = N;
  int C_INCX = INCX;
  int C_INCY = INCY;
  assert(X != NULL);
  assert(Y != NULL);
  assert(C_INCX > 0);
  assert(C_INCY > 0);
  BLAS_NAME(dcopy)(INTEGER(C_N), X, INTEGER(C_INCX), Y, INTEGER(C_INCY));
}

/** DGEMV - one of the matrix-vector operations y := alpha*A*x +
 * beta*y, or y := alpha*A'*x + xbeta*y,
 */
static inline void  DGEMV(const char * TRANS, int M, int N, double ALPHA, const double* const A, int LDA, const double* const X, int INCX, double BETA, double* Y, int INCY)
{
  int C_M = M;
  int C_N = N;
  double C_ALPHA = ALPHA;
  int C_LDA = LDA;
  int C_INCX = INCX;
  double C_BETA = BETA;
  int C_INCY = INCY;
  assert(C_LDA > 0);
  assert(C_M > 0);
  assert(C_N > 0);
  BLAS_NAME(dgemv)WITH_ORDER(LA_ORDER, T_TRANS(TRANS), INTEGER(C_M), INTEGER(C_N), DOUBLE(C_ALPHA), A, INTEGER(C_LDA), X, INTEGER(C_INCX), DOUBLE(C_BETA), Y, INTEGER(C_INCY));
}

/*  DGEMM  performs one of the matrix-matrix operations
 *
 *     C := alpha*op( A )*op( B ) + beta*C,
 *
 *  where  op( X ) is one of
 *
 *     op( X ) = X   or   op( X ) = X',*/
static inline void DGEMM(const char* TRANSA, const char* TRANSB, int M, int N, int K,
                         double ALPHA, double* A, int LDA, double* B, int LDB, double BETA, double* C, int LDC)
{
  int C_M = M;
  int C_N = N;
  int C_K = K;
  double C_ALPHA = ALPHA;
  int C_LDA = LDA;
  int C_LDB = LDB;
  double C_BETA = BETA;
  int C_LDC = LDC;
  BLAS_NAME(dgemm)WITH_ORDER(LA_ORDER, T_TRANS(TRANSA), T_TRANS(TRANSB), INTEGER(C_M), INTEGER(C_N),
                             INTEGER(C_K), DOUBLE(C_ALPHA), A, INTEGER(C_LDA), B, INTEGER(C_LDB), DOUBLE(C_BETA), C, INTEGER(C_LDC));
}

/* DDOT - the dot product of two vectors
*/
static inline double DDOT(int N, const double* const DX, int INCX, const double* const DY, int INCY)
{
  int C_N = N;
  int C_INCX = INCX;
  int C_INCY = INCY;
  return BLAS_NAME(ddot)(INTEGER(C_N), DX, INTEGER(C_INCX), DY, INTEGER(C_INCY));
}

/* DAXPY - time a vector plus a vector
*/
static inline void DAXPY(int N, double DA, double* DX, int INCX, double* DY, int INCY)
{
  int C_N = N;
  double C_DA = DA;
  int C_INCX = INCX;
  int C_INCY = INCY;
  BLAS_NAME(daxpy)(INTEGER(C_N), DOUBLE(C_DA), DX, INTEGER(C_INCX), DY, INTEGER(C_INCY));
}


/* DSCAL - a vector by a constant
*/
static inline void DSCAL(int N, double DA, double* DX, int INCX)
{
  int C_N = N;
  double C_DA = DA;
  int C_INCX = INCX;
  BLAS_NAME(dscal)(INTEGER(C_N), DOUBLE(C_DA), DX, INTEGER(C_INCX));
}

/* DPOTRF - compute the Cholesky factorization of a real symmetric
   positive definite matrix A
   */
static inline void DPOTRF(const char* UPLO, int N, double* A, int LDA, int* INFO)
{
  int C_N = N;
  int C_LDA = LDA;
  LAPACK_4(LAPACK_NAME(dpotrf), T_UPLO(UPLO), INTEGER(C_N), A , INTEGER(C_LDA), INTEGERP(INFO));
}

/* DPOTRI - compute the inverse of a real symmetric positive definite
   matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
   computed by DPOTRF
   */
/*static inline void DPOTRI(const char* UPLO, int N, double* A, int LDA, int* INFO)
  {
  int C_N = N;
  int C_LDA = LDA;
  LAPACK_4(LAPACK_NAME(dpotri), T_UPLO(UPLO), INTEGER(C_N), A, INTEGER(C_LDA), INTEGERP(INFO));
  }
  */
/* DGESV - compute the solution to a real system of linear equations A * X = B,
*/
static inline void DGESV(int N, int NRHS, double* A, int LDA, int* IPIV, double* B, int LDB, int* INFO)
{
  int C_N = N;
  int C_NRHS = NRHS;
  int C_LDA = LDA;
  int C_LDB = LDB;
  LAPACK_7(LAPACK_NAME(dgesv), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA), INTEGERP(IPIV), B, INTEGER(C_LDB), INTEGERP(INFO));
}


/* DGELS -
 *  DGELS solves overdetermined or underdetermined real linear systems
 *  involving an M-by-N matrix A, or its transpose, using a QR or LQ
 *  factorization of A.  It is assumed that A has full rank.
 */

#ifdef COMPLETE_LAPACK_LIBRARIES
#include "blaslapack.h"
static inline void DGELS(int M, int N, int NRHS, double* A, int LDA, double* B, int LDB, double* WORK, int LWORK, int* INFO)
{
  int C_M = M;
  int C_N = N;
  double * C_WORK = WORK;
  int C_NRHS = NRHS;
  int C_LDA = LDA;
  int C_LDB = LDB;
  int C_LWORK = LWORK;
  F77NAME(dgels)("N", FCAST(integer, C_M) , FCAST(integer, C_N), FCAST(integer, C_NRHS), FCASTP(double, A), FCAST(integer, C_LDA), FCASTP(double, B), FCAST(integer, C_LDB), FCASTP(double, C_WORK) , FCAST(integer, C_LWORK), FCASTP(integer, INFO));
}
#else
#include <stdio.h>
static inline void DGELS(int M, int N, int NRHS, double* A, int LDA, double* B, int LDB, double* WORK, int LWORK, int *INFO)
{
  *INFO = 1;
  fprintf(stderr, "DGELS not found. Please check you LAPACK/ATLAS installation.\n");
  fprintf(stderr, "%p", (void *)WORK);
  exit(EXIT_FAILURE);
}
#endif

/* DGESVD -
 * DGESVD  computes the singular value decomposition (SVD) of a real
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors.
 */
#ifdef COMPLETE_LAPACK_LIBRARIES
#include "blaslapack.h"
#define DGESVD(JOBU,JOBVT,M, N, A, LDA, S, U, LDU, VT, LDVT,WORK, LWORK, INFO ) \
  ({int C_M = M; \
   int C_N = N; \
   int C_LDA = LDA; \
   int C_LDU = LDU; \
   int C_LDVT = LDVT;\
   int C_LWORK = LWORK; \
   double * C_WORK = WORK; \
   F77NAME(dgesvd)(JOBU,JOBVT, FCAST(integer,C_M) , FCAST(integer,C_N),  FCASTP(double, A), FCAST(integer,C_LDA),FCASTP(double,S), \
     FCASTP(double,U),FCAST(integer, C_LDU), FCASTP(double,VT),FCAST(integer, C_LDVT),\
     FCASTP(double,C_WORK) ,  FCAST(integer,C_LWORK), FCAST(integer,INFO));\
   })
/* DGESVD -
 * DGESVD  computes the singular value decomposition (SVD) of a real
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors.
 */
#else
#include <stdio.h>
#define DGESVD(JOBU,JOBVT,M, N, A, LDA, S, U, LDU, VT, LDVT,WORK, LWORK, INFO ) \
  ({ fprintf(stderr, "dgesvd not found\n"); })
#endif


/* DTRTRS - solve a triangular system of the form  A * X = B or A**T * X = B,
*/
static inline void DTRTRS(const char* UPLO, const char* TRANS, const char* DIAG, int N, int NRHS, double* A, int LDA, double* B, int LDB, int* INFO)
{
  int C_N = N;
  int C_NRHS = NRHS;
  int C_LDA = LDA;
  int C_LDB = LDB;
  LAPACK_9_SIDED(LAPACK_NAME(dtrtrs), T_UPLO(UPLO), T_TRANS(TRANS), T_DIAG(DIAG), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA), B, INTEGER(C_LDB), INTEGERP(INFO));
}

/* DGETRF - LU factorization
*/
static inline void DGETRF(int M, int N, double* A, int LDA, int* IPIV, int* INFO)
{
  int C_M = M;
  int C_N = N;
  int C_LDA = LDA;
  LAPACK_5(LAPACK_NAME(dgetrf), INTEGER(C_M), INTEGER(C_N), A, INTEGER(C_LDA), INTEGERP(IPIV), INTEGERP(INFO));
}

/* DGETRS solves a system of linear equations
 *     A * X = B  or  A' * X = B
 *  with a general N-by-N matrix A using the LU factorization computed
 *  by DGETRF.
 */
static inline void DGETRS(const char* TRANS, int N, int NRHS, double* A, int LDA, int* IPIV, double* B, int LDB, int* INFO)
{
  int C_N = N;
  int C_NRHS = NRHS;
  int C_LDA = LDA;
  int C_LDB = LDB;
  LAPACK_8(LAPACK_NAME(dgetrs), T_TRANS(TRANS), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA), INTEGERP(IPIV), B, INTEGER(C_LDB), INTEGERP(INFO));
}

/* DGETRI - matrix inversion
*/
static inline void DGETRI(int N, double* A, int LDA, int* IPIV, int* INFO)
{
  int C_N = N;
  int C_LDA = LDA;
  printf("C_LDA = %i\n", C_LDA);
  LAPACK_4_W(LAPACK_NAME(dgetri), INTEGER(C_N), A, INTEGER(C_LDA), INTEGERP(IPIV), INTEGERP(INFO));
}


#endif /* USE_MKL */
#endif /* LA_H */
