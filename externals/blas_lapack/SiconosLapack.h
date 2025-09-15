/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef SiconosLAPACK_H
#define SiconosLAPACK_H

// Definition of the interface to cblas library.
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "SiconosBlas.h"
#include "SiconosConfig.h"

#if defined(__cplusplus)
extern "C" {
#endif

// Lapack implementation may differs in :
// - the name of the routines
// - the arguments of the routines
// - the presence of some routines (but only atlas seems to be lacking some functions and it's
// now forbidden in Siconos). Main differences for the arguments are :
//  - info parameter may be the last argument (apple framework) or a return value (lapacke,
//  mkl),
//  - a new argument is required to set the row or column major order for arrays (lapacke,
//  first arg.).
//  - some arg may be char (lapacke), char* (apple).
//
// At the time, we check for the following implementations (the check is done during cmake
// call): A - OpenBlas B - MKL intel C - Apple Accelerate framework D - lapacke from netlib
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

#elif defined(HAS_CLAPACK)
#include "SiconosAppleLapack.h"

// --- Intel MKL, OpenBlas and netlib lapacke ---
#else  //  defined(HAS_LAPACKE)
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
static inline void DGESV(lapack_int N, lapack_int NRHS, double* A, lapack_int LDA,
                         lapack_int* IPIV, double* B, lapack_int LDB, lapack_int* INFO) {
  lapack_int C_N = N;
  lapack_int C_NRHS = NRHS;
  lapack_int C_LDA = LDA;
  lapack_int C_LDB = LDB;
  lapack_int C_INFO = 0;
  WRAP_DGESV(LAPACK_NAME(dgesv), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA),
             INTEGERP(IPIV), B, INTEGER(C_LDB), INTEGER(C_INFO));
  *INFO = C_INFO;
}

static inline void DPOSV(char UPLO, lapack_int N, lapack_int NRHS, double* A, lapack_int LDA,
                         double* B, lapack_int LDB, lapack_int* INFO) {
  lapack_int C_N = N;
  lapack_int C_NRHS = NRHS;
  lapack_int C_LDA = LDA;
  lapack_int C_LDB = LDB;
  lapack_int C_INFO = 0;
  WRAP_DPOSV(LAPACK_NAME(dposv), CHAR(UPLO), INTEGER(C_N), INTEGER(C_NRHS), A, INTEGER(C_LDA),
             B, INTEGER(C_LDB), INTEGER(C_INFO));
  *INFO = C_INFO;
}

/* DGETRF - LU factorization */
static inline void DGETRF(lapack_int M, lapack_int N, double* A, lapack_int LDA,
                          lapack_int* IPIV, lapack_int* INFO) {
  lapack_int C_M = M;
  lapack_int C_N = N;
  lapack_int C_LDA = LDA;
  lapack_int C_INFO = 0;
  WRAP_DGETRF(LAPACK_NAME(dgetrf), INTEGER(C_M), INTEGER(C_N), A, INTEGER(C_LDA),
              INTEGERP(IPIV), INTEGER(C_INFO));
  *INFO = C_INFO;
}

/* DGETRI - matrix inversion
 */
static inline void DGETRI(lapack_int N, double* A, lapack_int LDA, lapack_int* IPIV,
                          lapack_int* INFO) {
  lapack_int C_N = N;
  lapack_int C_LDA = LDA;
  lapack_int C_INFO = 0;
  WRAP_DGETRI(LAPACK_NAME(dgetri), INTEGER(C_N), A, INTEGER(C_LDA), INTEGERP(IPIV),
              INTEGER(C_INFO));
  *INFO = C_INFO;
}

/* DSYTRF - LDLT factorization */
static inline void DSYTRF(char UPLO, lapack_int N, double* A, lapack_int LDA, lapack_int* IPIV,
                          lapack_int* INFO) {
  lapack_int C_N = N;
  lapack_int C_LDA = LDA;
  lapack_int C_INFO = 0;
  WRAP_DSYTRF(LAPACK_NAME(dsytrf), CHAR(UPLO), INTEGER(C_N), A, INTEGER(C_LDA), INTEGERP(IPIV),
              INTEGER(C_INFO));
  *INFO = C_INFO;
}

/*
  DGESVD - The routine computes the singular value decomposition (SVD) of a rectangular real
  matrix A, optionally the left and/or right singular vectors. Note Franck : there is no common
  interface for the different lapack implementations. Either work/lwork must be provided (as in
  F77 interf, for internal memory management) or not. Some interfaces (lapacke from netlib,
  OpenBlas) provide both (middle and high level interfaces, see
  http://www.netlib.org/lapack/lapacke.html). Actually only Apple framework is missing the "no
  work/lwork" interface. So, two different cases :
  - fortran like interface, that needs a properly allocated work.
  - lapacke high-level interface, with work == superb, memory allocation hide from used.
*/
static inline void DGESVD(char JOBU, char JOBVT, lapack_int M, lapack_int N, double* A,
                          lapack_int LDA, double* S, double* U, lapack_int LDU, double* VT,
                          lapack_int LDVT, double* superb, lapack_int* INFO) {
  lapack_int C_M = M;
  lapack_int C_N = N;
  lapack_int C_LDA = LDA;
  lapack_int C_LDU = LDU;
  lapack_int C_LDVT = LDVT;
  lapack_int C_INFO = 0;
  WRAP_DGESVD(LAPACK_NAME(dgesvd), CHAR(JOBU), CHAR(JOBVT), INTEGER(C_M), INTEGER(C_N), A,
              INTEGER(C_LDA), S, U, INTEGER(C_LDU), VT, INTEGER(C_LDVT), superb,
              INTEGER(C_INFO));
  *INFO = C_INFO;
}

/* DGETRS solves a system of linear equations
 *     A * X = B  or  A' * X = B
 *  with a general N-by-N matrix A using the LU factorization computed
 *  by DGETRF.
 */
static inline void DGETRS(char TRANS, lapack_int N, lapack_int NRHS, double* A, lapack_int LDA,
                          lapack_int* IPIV, double* B, lapack_int LDB, lapack_int* INFO) {
  lapack_int C_N = N;
  lapack_int C_NRHS = NRHS;
  lapack_int C_LDA = LDA;
  lapack_int C_LDB = LDB;
  lapack_int C_INFO = 0;
  WRAP_DGETRS(LAPACK_NAME(dgetrs), CHAR(TRANS), INTEGER(C_N), INTEGER(C_NRHS), A,
              INTEGER(C_LDA), INTEGERP(IPIV), B, INTEGER(C_LDB), INTEGER(C_INFO));
  *INFO = C_INFO;
}

/*
  DGELS - The routine solves overdetermined or underdetermined real linear
  systems involving a rectangular matrix A, or its transpose,
  using a QR or LQ factorization of A.
  The routine assumes that A has full rank.
*/
static inline void DGELS(char trans, lapack_int M, lapack_int N, lapack_int NRHS, double* A,
                         lapack_int LDA, double* B, lapack_int LDB, lapack_int* INFO) {
  lapack_int C_M = M;
  lapack_int C_N = N;
  lapack_int C_NRHS = NRHS;
  lapack_int C_LDA = LDA;
  lapack_int C_LDB = LDB;
  lapack_int C_INFO = 0;
  WRAP_DGELS(LAPACK_NAME(dgels), CHAR(trans), INTEGER(C_M), INTEGER(C_N), INTEGER(C_NRHS), A,
             INTEGER(C_LDA), B, INTEGER(C_LDB), INTEGER(C_INFO));
  *INFO = C_INFO;
}

/* DPOTRF - compute the Cholesky factorization of a real symmetric
   positive definite matrix A
*/
static inline void DPOTRF(char UPLO, lapack_int N, double* A, lapack_int LDA,
                          lapack_int* INFO) {
  lapack_int C_N = N;
  lapack_int C_LDA = LDA;
  lapack_int C_INFO = 0;
  WRAP_DPOTRF(LAPACK_NAME(dpotrf), CHAR(UPLO), INTEGER(C_N), A, INTEGER(C_LDA),
              INTEGER(C_INFO));
  *INFO = C_INFO;
}
/* DPOTRS solves a system of linear equations
 *     A * X = B
 *  with a symmetrix definite positive N-by-N matrix A using the Cholesky factorization
 * computed by DPOTRF.
 */
static inline void DPOTRS(char UPLO, lapack_int N, lapack_int NRHS, double* A, lapack_int LDA,
                          double* B, lapack_int LDB, lapack_int* INFO) {
  lapack_int C_N = N;
  lapack_int C_NRHS = NRHS;
  lapack_int C_LDA = LDA;
  lapack_int C_LDB = LDB;
  lapack_int C_INFO = 0;
  WRAP_DPOTRS(LAPACK_NAME(dpotrs), CHAR(UPLO), INTEGER(C_N), INTEGER(C_NRHS), A,
              INTEGER(C_LDA), B, INTEGER(C_LDB), INTEGER(C_INFO));
  *INFO = C_INFO;
}
/* DSYTRS solves a system of linear equations
 *     A * X = B
 *  with a general symmetric N-by-N matrix A using the LDLT factorization computed
 *  by DSYTRF.
 */
static inline void DSYTRS(char UPLO, lapack_int N, lapack_int NRHS, double* A, lapack_int LDA,
                          lapack_int* IPIV, double* B, lapack_int LDB, lapack_int* INFO) {
  lapack_int C_N = N;
  lapack_int C_NRHS = NRHS;
  lapack_int C_LDA = LDA;
  lapack_int C_LDB = LDB;
  lapack_int C_INFO = 0;
  WRAP_DSYTRS(LAPACK_NAME(dsytrs), CHAR(UPLO), INTEGER(C_N), INTEGER(C_NRHS), A,
              INTEGER(C_LDA), INTEGERP(IPIV), B, INTEGER(C_LDB), INTEGER(C_INFO));
  *INFO = C_INFO;
}
/* DTRTRS - solve a triangular system of the form  A * X = B or A**T * X = B,
 */
static inline void DTRTRS(char UPLO, char TRANS, char DIAG, lapack_int N, lapack_int NRHS,
                          double* A, lapack_int LDA, double* B, lapack_int LDB,
                          lapack_int* INFO) {
  lapack_int C_N = N;
  lapack_int C_NRHS = NRHS;
  lapack_int C_LDA = LDA;
  lapack_int C_LDB = LDB;
  lapack_int C_INFO = 0;
  WRAP_DTRTRS(LAPACK_NAME(dtrtrs), CHAR(UPLO), CHAR(TRANS), CHAR(DIAG), INTEGER(C_N),
              INTEGER(C_NRHS), A, INTEGER(C_LDA), B, INTEGER(C_LDB), INTEGER(C_INFO));
  *INFO = C_INFO;
}

#if defined(__cplusplus)
}
#endif

#endif  // SICONOSLAPACK_H
