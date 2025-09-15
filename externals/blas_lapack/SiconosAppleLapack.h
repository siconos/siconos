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

#ifndef SiconosAppleLAPACK_H
#define SiconosAppleLAPACK_H

// IWYU pragma: private, include "SiconosLapack.h"

// #include "SiconosBlas.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define FCAST(T, X) (T*)(&X)
#define FCASTP(T, X) (T*)X

#include <clapack.h>
#define LAPACK_NAME(N) N##_

#define LA_TRANS 'T'
#define LA_NOTRANS 'N'
#define LA_UP 'U'
#define LA_LO 'L'
#define LA_NONUNIT 'N'
#define LA_UNIT 'U'
#define INTEGER(X) FCAST(__CLPK_integer, X)
#define INTEGERP(X) FCASTP(__CLPK_integer, X)
#define CHAR(X) FCAST(char, X)

// --- DGESVD ---
// Note FP : we need to call WRAP_DGESVD two times, one to find optimal lwork
// and one for the real computation.
// We may add a test on info value to call lwork only once?
#define WRAP_DGESVD(F, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, INFO) \
  __CLPK_integer lwork = -1;                                                    \
  __CLPK_doublereal* work = (__CLPK_doublereal*)malloc(sizeof(*work));          \
  F(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, work, &lwork, INFO);          \
  lwork = (__CLPK_integer)work[0];                                              \
  work = (__CLPK_doublereal*)realloc(work, (size_t)lwork * sizeof(*work));      \
  F(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, work, &lwork, INFO);          \
  free(work);

// --- DGETRS ---
#define WRAP_DGETRS(F, A1, A2, A3, A4, A5, A6, A7, A8, INFO) \
  F(A1, A2, A3, A4, A5, A6, A7, A8, INFO)

// --- DPOTRS ---
#define WRAP_DPOTRS(F, A1, A2, A3, A4, A5, A6, A7, INFO) F(A1, A2, A3, A4, A5, A6, A7, INFO)

// --- DSYTRS ---
#define WRAP_DSYTRS(F, A1, A2, A3, A4, A5, A6, A7, A8, INFO) \
  F(A1, A2, A3, A4, A5, A6, A7, A8, INFO)

// --- DGESV ---
#define WRAP_DGESV(F, A1, A2, A3, A4, A5, A6, A7, INFO) F(A1, A2, A3, A4, A5, A6, A7, INFO)

// --- DPOSV ---
#define WRAP_DPOSV(F, A1, A2, A3, A4, A5, A6, A7, INFO) F(A1, A2, A3, A4, A5, A6, A7, INFO)

// --- DGELS ---
#define WRAP_DGELS(F, A1, A2, A3, A4, A5, A6, A7, A8, INFO)                \
  __CLPK_integer lwork = -1;                                               \
  __CLPK_doublereal* work = (__CLPK_doublereal*)malloc(sizeof(*work));     \
  F(A1, A2, A3, A4, A5, A6, A7, A8, work, &lwork, INFO);                   \
  lwork = (__CLPK_integer)work[0];                                         \
  work = (__CLPK_doublereal*)realloc(work, (size_t)lwork * sizeof(*work)); \
  F(A1, A2, A3, A4, A5, A6, A7, A8, work, &lwork, INFO);                   \
  free(work);

// --- DGETRI ---
#define WRAP_DGETRI(F, A1, A2, A3, A4, INFO)                                    \
  __CLPK_integer lwork = -1;                                                    \
  __CLPK_doublereal* C_WORK;                                                    \
  C_WORK = (__CLPK_doublereal*)malloc(sizeof *C_WORK);                          \
  assert(C_WORK);                                                               \
  F(A1, A2, A3, A4, C_WORK, &lwork, INFO);                                      \
  lwork = (__CLPK_integer)(C_WORK[0]);                                          \
  C_WORK = (__CLPK_doublereal*)realloc(C_WORK, (size_t)lwork * sizeof *C_WORK); \
  F(A1, A2, A3, A4, C_WORK, &lwork, INFO);                                      \
  free(C_WORK);

// --- DGETRF ---
#define WRAP_DGETRF(F, A1, A2, A3, A4, A5, INFO) F(A1, A2, A3, A4, A5, INFO)

// --- DPOTRF ---
#define WRAP_DPOTRF(F, A1, A2, A3, A4, INFO) F(A1, A2, A3, A4, INFO);

// --- DSYTRF ---
#define WRAP_DSYTRF(F, A1, A2, A3, A4, A5, INFO)                     \
  lapack_int lwork = -1;                                             \
  double* C_WORK;                                                    \
  C_WORK = (double*)malloc(sizeof *C_WORK);                          \
  assert(C_WORK);                                                    \
  F(A1, A2, A3, A4, A5, C_WORK, &lwork, INFO);                       \
  lwork = (lapack_int)(C_WORK[0]);                                   \
  C_WORK = (double*)realloc(C_WORK, (size_t)lwork * sizeof *C_WORK); \
  F(A1, A2, A3, A4, A5, C_WORK, &lwork, INFO);                       \
  free(C_WORK);

// --- DTRTRS ---
#define WRAP_DTRTRS(F, A1, A2, A3, A4, A5, A6, A7, A8, A9, INFO) \
  F(A1, A2, A3, A4, A5, A6, A7, A8, A9, INFO)

#endif  // SICONOSAPPLELAPACK_H
