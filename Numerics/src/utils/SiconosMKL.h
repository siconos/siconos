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

/* This files describes interfaces for blas/lapack functions using the MKL intel library.
   It must be used if Numerics is compiled/linked using "USE_MKL=TRUE" in cmake process.
   It's used in place of atlas or similar libraries.

   Note Franck: mkl is supposed to work either with intel or gnu compilers.
*/

// About blas, lapack in mkl see
// http://software.intel.com/sites/products/documentation/hpc/mkl/mklman/index.htm
// http://software.intel.com/en-us/articles/intel-math-kernel-library-intel-mkl-blas-cblas-and-lapack-compilinglinking-functions-fortran-and-cc-calls/
//
//

#ifndef SiconosMK_H
#define SiconosMKL

#define min(a,b) ((a)>(b)?(b):(a))

// ================ BLAS FUNCTIONS ================
// Those functions are declared in mkl_cblas.h
// and use prefix cblas_
#include <mkl_cblas.h>
#define BLAS_NAME(N) cblas_##N
#ifndef LA_ORDER
#define LA_ORDER CblasColMajor
#endif
#define LA_TRANS CblasTrans
#define LA_NOTRANS CblasNoTrans

/* DNRM2 - the euclidean norm of a vector via the function name, so
   that DNRM2 := sqrt( x'*x )
*/
#define DNRM2(N, X, INCX)                                \
  ({  assert(N>0);                                       \
      assert(INCX>0);                                    \
      assert(X!=NULL);                                     \
      BLAS_NAME(dnrm2)(N, X, INCX);  \
  })

/* DCOPY - a vector, x, to a vector, y
*/
#define DCOPY(N, X, INCX, Y, INCY)                                      \
  ({ assert(X!=NULL);                                                    \
     assert(Y!=NULL);                                                    \
     assert(INCX>0);                                                   \
     assert(INCY>0);                                                   \
     BLAS_NAME(dcopy)(N, X, INCX, Y, INCY); \
  })

/** DGEMV - one of the matrix-vector operations y := alpha*A*x +
   beta*y, or y := alpha*A'*x + xbeta*y,
*/
#define DGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)       \
  ({ assert (LDA > 0);                                                 \
     assert (M > 0);                                                   \
     assert (N > 0);                                                   \
     BLAS_NAME(dgemv)(LA_ORDER, TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY); \
  })

/* DDOT - the dot product of two vectors */
#define DDOT(N, DX, INCX, DY, INCY) \
  ({ BLAS_NAME(ddot)(N, DX, INCX, DY, INCY);  })

/*DAXPY - time a vector plus a vector */
#define DAXPY(N, DA, DX, INCX, DY, INCY) \
  ({ BLAS_NAME(daxpy)(N, DA,DX,INCX,DY,INCY); })

/* DSCAL - a vector by a constant */
#define DSCAL(N,DA,DX,INCX) \
  ({ BLAS_NAME(dscal)(N, DA, DX, INCX); })

/*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',*/
#define DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) \
  ({ BLAS_NAME(dgemm)(LA_ORDER,TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC); })

//  ================ LAPACK FUNCTIONS ================
// Those functions are declared in mkl_lapacke.h
// and use prefix LAPACKE_
#include <mkl_lapacke.h>
#define LAPACK_NAME(N) LAPACKE_##N
// Note Franck : We need to redefine the equivalent for CBlasTrans and NoTrans since in clapack interface of MKL, they are declared as char, not int ...
#define CLA_TRANS 'T'
#define CLA_NOTRANS 'N'
#define LA_UP 'U'
#define LA_LO 'L'
#define LA_NONUNIT 'N' // 'N' or CblasNonUnit
#define LA_UNIT 'U'    // 'U' or CblasUnit


/* DPOTRF - compute the Cholesky factorization of a real symmetric
   positive definite matrix A
 */
#define DPOTRF(UPLO, N, A, LDA, INFO) \
  ({ INFO = LAPACK_NAME(dpotrf)(LA_ORDER, UPLO, N,A,LDA); })

/* DTRTRS - solve a triangular system of the form  A * X = B or A**T * X = B, */
#define DTRTRS(UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO ) \
  ({ INFO = LAPACK_NAME(dtrtrs)(LA_ORDER, UPLO,TRANS,DIAG,N,NRHS,A,LDA,B,LDB); })

/* DGETRF - LU factorization */
#define DGETRF(M,N,A,LDA,IPIV,INFO) \
  ({ INFO = LAPACK_NAME(dgetrf)(LA_ORDER,M,N,A,LDA,IPIV); })

/* DGESV - compute the solution to a real system of linear equations A * X = B, */
#define DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO ) \
  ({ INFO = LAPACK_NAME(dgesv)(LA_ORDER,N, NRHS, A, LDA, IPIV, B, LDB);  })

/* DGELS -
 *  DGELS solves overdetermined or underdetermined real linear systems
 *  involving an M-by-N matrix A, or its transpose, using a QR or LQ
 *  factorization of A.  It is assumed that A has full rank.
 */
#define DGELS( M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO  ) \
  ({ INFO = LAPACK_NAME(dgels)(LA_ORDER,'N',M,N,NRHS,A,LDA,B,LDB);  })

/* DGETRI - matrix inversion */
#define DGETRI(N,A,LDA,IPIV,INFO) \
  ({ INFO = LAPACK_NAME(dgetri)(LA_ORDER,N,A,LDA,IPIV); })

/* DGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
 */
#define DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO) \
    ({ INFO = LAPACK_NAME(dgetrs)(LA_ORDER,TRANS,N,NRHS,A,LDA,IPIV,B,LDB); })

/* DGESVD -
 * DGESVD  computes the singular value decomposition (SVD) of a real
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors.
 */
#define DGESVD(JOBU,JOBVT,M, N, A, LDA, S, U, LDU, VT, LDVT,WORK, LWORK, INFO ) \
  ({ INFO = LAPACK_NAME(dgesvd)(LA_ORDER,JOBU,JOBVT,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK); })

#endif



