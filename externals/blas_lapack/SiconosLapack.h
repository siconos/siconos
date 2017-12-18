/* Siconos-Numerics, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#ifndef SiconosLAPACK_H
#define SiconosLAPACK_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

// Definition of the interface to cblas library. 
#include "SiconosBlas.h"

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#if defined(__cplusplus)
extern "C"
{
#endif


// Lapack implementation may differs in : 
// - the name of the routines 
// - the arguments of the routines
// - the presence of some routines (only atlas seems to be lacking some functions).
// Main differences for the arguments are :
//  - info parameter may be the last argument (apple framework) or a return value (lapacke, mkl, atlas),
//  - a new argument is required to set the row or column major order for arrays (lapacke, atlas, first arg.).
//  - some arg may be char (lapacke), char* (apple) or enum (atlas). 
//
// At the time, we check for the following implementations (the check is done during cmake call):
// A - OpenBlas
// B - MKL intel
// C - Atlas
// D - Apple Accelerate framework
// E - lapacke from netlib
//
// A, B and E have common interf. (lapacke)
// D is more or less close to clapack.c
// C looks like lapacke but is not complete. 
// 
// Thus, we define some macros in the headers below, to set a common interface for Siconos
// calls to lapack functions. 
// The HAS_... variables are set by cmake and defined in SiconosConfig.h and
// the interfaces to SiconosLapack functions are listed at the bottom of this file. 
//
// Each include defines : 
// - Headers and routines naming conventions for the different Lapack implementations 
// - WRAP_XXX functions
//


// ---  Apple Accelerate framework ---
#if defined(HAS_ACCELERATE)
#include "SiconosAppleLapack.h"

 // --- MATLAB ---
#elif defined(HAS_MATLAB_LAPACK)
#include "SiconosMATLABLapack.h"

// --- Atlas ---
#elif defined(HAS_ATLAS_LAPACK) || !defined(HAS_LAPACKE)
#include "SiconosAtlasLapack.h"

#elif defined(HAS_CLAPACK)
#include "SiconosAppleLapack.h"

// --- Intel MKL, OpenBlas and netlib lapacke --- 
#else //  defined(HAS_LAPACKE) 
#include "SiconosLapacke.h"
#endif

#ifndef lapack_int
#define lapack_int int
#endif

/*
  DGESV - The routine solves the system of linear equations for X:
  A*X = B
  where
  A is a square matrix.
  The columns of matrix B are individual right-hand sides.
  The columns of X are the corresponding solutions.
*/
  static inline void DGESV(lapack_int N, lapack_int NRHS, double* A, lapack_int LDA, lapack_int* IPIV, double* B, lapack_int LDB, lapack_int* INFO)
  {
    lapack_int C_N = N;
    lapack_int C_NRHS = NRHS;
    lapack_int C_LDA = LDA;
    lapack_int C_LDB = LDB;
    lapack_int C_INFO = 0;
    WRAP_DGESV(LAPACK_NAME(dgesv), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA), INTEGERP(IPIV), B, INTEGER(C_LDB), INTEGER(C_INFO));
    *INFO = C_INFO;
  }

/* DGETRF - LU factorization */
  static inline void DGETRF(lapack_int M, lapack_int N, double* A, lapack_int LDA, lapack_int* IPIV, lapack_int* INFO)
  {
    lapack_int C_M = M;
    lapack_int C_N = N;
    lapack_int C_LDA = LDA;
    lapack_int C_INFO = 0;
    WRAP_DGETRF(LAPACK_NAME(dgetrf), INTEGER(C_M), INTEGER(C_N), A, INTEGER(C_LDA), INTEGERP(IPIV), INTEGER(C_INFO));
    *INFO = C_INFO;
  }

/* DGETRI - matrix inversion
 */
  static inline void DGETRI(lapack_int N, double* A, lapack_int LDA, lapack_int* IPIV, lapack_int* INFO)
  {
    lapack_int C_N = N;
    lapack_int C_LDA = LDA;
    lapack_int C_INFO = 0;
    WRAP_DGETRI(LAPACK_NAME(dgetri), INTEGER(C_N), A, INTEGER(C_LDA), INTEGERP(IPIV), INTEGER(C_INFO));
    *INFO = C_INFO;
  }


#if defined(HAS_ATLAS_LAPACK)

  static inline void DGESVD(char JOBU, char JOBVT, lapack_int M, lapack_int N, double * A, lapack_int LDA, double * S, double * U, lapack_int LDU, double * VT, lapack_int LDVT, double * superb, lapack_int* INFO)
  {
    /* int C_M = M; */
    /* int C_N = N; */
    /* int C_LDA = LDA; */
    /* int C_LDU = LDU; */
    /* int C_LDVT = LDVT; */
    /* int C_INFO = 0; */
    WRAP_DGESVD(LAPACK_NAME(dgesvd), CHAR(JOBU), CHAR(JOBVT), INTEGER(C_M), INTEGER(C_N), A, INTEGER(C_LDA), S, U, INTEGER(C_LDU), VT, INTEGER(C_LDVT), superb, INTEGER(C_INFO));
//    *INFO = C_INFO;
  }

/* DGETRS solves a system of linear equations
 *     A * X = B  or  A' * X = B
 *  with a general N-by-N matrix A using the LU factorization computed
 *  by DGETRF.
 */
  static inline void DGETRS(const enum CBLAS_TRANSPOSE TRANS, lapack_int N, lapack_int NRHS, double* A, lapack_int LDA, lapack_int* IPIV, double* B, lapack_int LDB, lapack_int* INFO)
  {
    lapack_int C_N = N;
    lapack_int C_NRHS = NRHS;
    lapack_int C_LDA = LDA;
    lapack_int C_LDB = LDB;
    lapack_int C_INFO = 0;
    WRAP_DGETRS(LAPACK_NAME(dgetrs), TRANS, INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA), INTEGERP(IPIV), B, INTEGER(C_LDB), INTEGER(C_INFO));
    *INFO = C_INFO;
  }

/*
  DGELS - The routine solves overdetermined or underdetermined real linear
  systems involving a rectangular matrix A, or its transpose, 
  using a QR or LQ factorization of A. 
  The routine assumes that A has full rank.
*/
  static inline void DGELS(const enum CBLAS_TRANSPOSE trans, lapack_int M, lapack_int N, lapack_int NRHS, double* A, lapack_int LDA, double* B, lapack_int LDB, lapack_int* INFO)
  {
    /* lapack_int C_M = M; */
    /* lapack_int C_N = N; */
    /* lapack_int C_NRHS = NRHS; */
    /* lapack_int C_LDA = LDA; */
    /* lapack_int C_LDB = LDB; */
    /* lapack_int C_INFO = 0 =*INFO; */
    WRAP_DGELS(LAPACK_NAME(dgels),trans, INTEGER(C_M), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA), B, INTEGER(C_LDB),INTEGER(C_INFO));
    /* *INFO = C_INFO; */
  }

  /* DPOTRF - compute the Cholesky factorization of a real symmetric
     positive definite matrix A
  */
  static inline void DPOTRF(const enum ATLAS_UPLO UPLO, lapack_int N, double* A, lapack_int LDA, lapack_int* INFO)
  {
    lapack_int C_N = N;
    lapack_int C_LDA = LDA;
    lapack_int C_INFO = 0;
    WRAP_DPOTRF(LAPACK_NAME(dpotrf), UPLO, INTEGER(C_N), A , INTEGER(C_LDA), INTEGER(C_INFO));
    *INFO = C_INFO;
  }

  /* DTRTRS - solve a triangular system of the form  A * X = B or A**T * X = B,
   */
  static inline void DTRTRS(const enum ATLAS_UPLO UPLO, const enum CBLAS_TRANSPOSE TRANS, const enum CBLAS_DIAG  DIAG, lapack_int N, lapack_int NRHS, double* A, lapack_int LDA, double* B, lapack_int LDB, lapack_int* INFO)
  {
    *INFO = clapack_dtrtrs(CblasColMajor, CblasLeft, UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB);
//  WRAP_DTRTRS(LAPACK_NAME(dtrtrs), CHAR(UPLO), CHAR(TRANS), CHAR(DIAG), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA), B, INTEGER(C_LDB), INTEGER(C_INFO));
  }


#else // No atlas

  
/*
  DGESVD - The routine computes the singular value decomposition (SVD) of a rectangular real matrix A, optionally the left and/or right singular vectors.
  Note Franck : there is no common interface for the different lapack implementations.
  Either work/lwork must be provided (as in F77 interf, for internal memory management) or not.
  Some interfaces (lapacke from netlib, OpenBlas) provide both (middle and high level interfaces, see http://www.netlib.org/lapack/lapacke.html).
  Actually only Apple framework is missing the "no work/lwork" interface, while atlas has no dgesvd at all.
  So, two different cases :
  - fortran like interface, that needs a properly allocated work.
  - lapacke high-level interface, with work == superb, memory allocation hide from used.
*/
  static inline void DGESVD(char JOBU, char JOBVT, lapack_int M, lapack_int N, double * A, lapack_int LDA, double * S, double * U, lapack_int LDU, double * VT, lapack_int LDVT, double * superb, lapack_int* INFO)
  {
    lapack_int C_M = M;
    lapack_int C_N = N;
    lapack_int C_LDA = LDA;
    lapack_int C_LDU = LDU;
    lapack_int C_LDVT = LDVT;
    lapack_int C_INFO = 0;
    WRAP_DGESVD(LAPACK_NAME(dgesvd), CHAR(JOBU), CHAR(JOBVT), INTEGER(C_M), INTEGER(C_N), A, INTEGER(C_LDA), S, U, INTEGER(C_LDU), VT, INTEGER(C_LDVT), superb, INTEGER(C_INFO));
    *INFO = C_INFO;
  }

/* DGETRS solves a system of linear equations
 *     A * X = B  or  A' * X = B
 *  with a general N-by-N matrix A using the LU factorization computed
 *  by DGETRF.
 */
  static inline void DGETRS(char TRANS, lapack_int N, lapack_int NRHS, double* A, lapack_int LDA, lapack_int* IPIV, double* B, lapack_int LDB, lapack_int* INFO)
  {
    lapack_int C_N = N;
    lapack_int C_NRHS = NRHS;
    lapack_int C_LDA = LDA;
    lapack_int C_LDB = LDB;
    lapack_int C_INFO = 0;
    WRAP_DGETRS(LAPACK_NAME(dgetrs), CHAR(TRANS), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA), INTEGERP(IPIV), B, INTEGER(C_LDB), INTEGER(C_INFO));
    *INFO = C_INFO;
  }


/*
  DGELS - The routine solves overdetermined or underdetermined real linear
  systems involving a rectangular matrix A, or its transpose, 
  using a QR or LQ factorization of A. 
  The routine assumes that A has full rank.
*/
  static inline void DGELS(char trans, lapack_int M, lapack_int N, lapack_int NRHS, double* A, lapack_int LDA, double* B, lapack_int LDB, lapack_int* INFO)
  {
    lapack_int C_M = M;
    lapack_int C_N = N;
    lapack_int C_NRHS = NRHS;
    lapack_int C_LDA = LDA;
    lapack_int C_LDB = LDB;
    lapack_int C_INFO = 0;
    WRAP_DGELS(LAPACK_NAME(dgels),CHAR(trans), INTEGER(C_M), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA), B, INTEGER(C_LDB),INTEGER(C_INFO));
    *INFO = C_INFO;
  }

  /* DPOTRF - compute the Cholesky factorization of a real symmetric
     positive definite matrix A
  */
  static inline void DPOTRF(char UPLO, lapack_int N, double* A, lapack_int LDA, lapack_int* INFO)
  {
    lapack_int C_N = N;
    lapack_int C_LDA = LDA;
    lapack_int C_INFO = 0;
    WRAP_DPOTRF(LAPACK_NAME(dpotrf), CHAR(UPLO), INTEGER(C_N), A , INTEGER(C_LDA), INTEGER(C_INFO));
    *INFO = C_INFO;
  }

  /* DTRTRS - solve a triangular system of the form  A * X = B or A**T * X = B,
   */
  static inline void DTRTRS(char UPLO, char TRANS, char DIAG, lapack_int N, lapack_int NRHS, double* A, lapack_int LDA, double* B, lapack_int LDB, lapack_int* INFO)
  {
    lapack_int C_N = N;
    lapack_int C_NRHS = NRHS;
    lapack_int C_LDA = LDA;
    lapack_int C_LDB = LDB;
    lapack_int C_INFO = 0;
    WRAP_DTRTRS(LAPACK_NAME(dtrtrs), CHAR(UPLO), CHAR(TRANS), CHAR(DIAG), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA), B, INTEGER(C_LDB), INTEGER(C_INFO));
    *INFO = C_INFO;
  }
#endif // No atlas


#if defined(__cplusplus)
}
#endif

#pragma GCC diagnostic pop

#endif // SICONOSLAPACK_H
