/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#ifndef SiconosLAPACKE_H
#define SiconosLAPACKE_H

// IWYU pragma: private, include "SiconosLapack.h"
//#include "SiconosBlas.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

// -------- Headers and routines naming conventions for the different Lapack implementations --------

// --- Intel MKL Header ---
#if defined(HAS_MKL_LAPACKE)
#include <mkl_lapacke.h>
#else
// Standard lapacke header
#include <lapacke.h>
#endif

// Name of the routines
#define LAPACK_NAME(N) LAPACKE_##N

#define LA_TRANS 'T'
#define LA_NOTRANS 'N'
#define LA_UP 'U'
#define LA_LO 'L'
#define LA_NONUNIT 'N'
#define LA_UNIT 'U'

#define INTEGER(X) X
#define INTEGERP(X) X
#define CHAR(X) X

#ifndef lapack_int
#define lapack_int int
#endif

// --- DGESVD ---
#if defined(HAS_LAPACK_dgesvd)
#define WRAP_DGESVD(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,INFO)      \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12)
#else
#define WRAP_DGESVD(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,INFO)      \
  fprintf(stderr, "Your lapack version misses dgesvd function.\n");
#endif

// --- DGETRS ---
#define WRAP_DGETRS(F,A1,A2,A3,A4,A5,A6,A7,A8,INFO)   \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5,A6,A7,A8)

// --- DPOTRS ---
#define WRAP_DPOTRS(F,A1,A2,A3,A4,A5,A6,A7,INFO) \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5,A6,A7)

// --- DSYTRS ---
#define WRAP_DSYTRS(F,A1,A2,A3,A4,A5,A6,A7,A8,INFO) \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5,A6,A7,A8)


// --- DGESV ---
#define WRAP_DGESV(F,A1,A2,A3,A4,A5,A6,A7,INFO)   \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5,A6,A7)

// --- DPOSV ---
#define WRAP_DPOSV(F,A1,A2,A3,A4,A5,A6,A7,INFO) \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5,A6,A7)

// --- DGELS ---
#if defined(HAS_LAPACK_dgels)
#define WRAP_DGELS(F,A1,A2,A3,A4,A5,A6,A7,A8,INFO)    \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5,A6,A7,A8)
#else
#define WRAP_DGELS(F,A1,A2,A3,A4,A5,A6,A7,A8,INFO)                      \
  fprintf(stderr, "Your lapack version misses dgels function.\n");
#endif

// --- DGETRI ---
#define WRAP_DGETRI(F,A1,A2,A3,A4,INFO)         \
  INFO = F(CblasColMajor,A1,A2,A3,A4)


// --- DGETRF ---
#define WRAP_DGETRF(F,A1,A2,A3,A4,A5,INFO)  \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5)

// --- DPOTRF ---
#define WRAP_DPOTRF(F,A1,A2,A3,A4,INFO)  \
  INFO = F(CblasColMajor,A1,A2,A3,A4)

// --- DSYTRF ---
#define WRAP_DSYTRF(F,A1,A2,A3,A4,A5,INFO)      \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5)

// --- DTRTRS ---
#if defined(HAS_LAPACK_dtrtrs)
#define WRAP_DTRTRS(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,INFO)  \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5,A6,A7,A8,A9)
#else
#define WRAP_DTRTRS(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,INFO)                \
  fprintf(stderr, "Your lapack version misses dtrtrs function.\n");
#endif

#endif // SICONOSLAPACKE_H
