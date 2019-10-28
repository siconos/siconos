/* Siconos-Numerics, Copyright INRIA 2005-2018.
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

#ifndef SiconosAtlasLAPACK_H
#define SiconosAtlasLAPACK_H

#include "SiconosBlas.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

// -------- Headers and routines naming conventions for the different Lapack implementations --------

#include <clapack.h>
#define LAPACK_NAME(N) clapack_##N
#define LA_TRANS CblasTrans
#define LA_NOTRANS CblasNoTrans
#define LA_UP CblasUpper
#define LA_LO CblasLower
#define LA_NONUNIT CblasNonUnit
#define LA_UNIT CblasUnit
#define INTEGER(X) X
#define INTEGERP(X) X
#define CHAR(X) X

int clapack_dtrtrs(const enum ATLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum ATLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans, const enum CBLAS_DIAG Diag, int n, int nrhs, double *a, int lda, double *b, int ldb);

// --- DGESVD ---
#define WRAP_DGESVD(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,INFO)      \
  fprintf(stderr, "Your lapack (atlas) version misses dgesvd function.\n");

// --- DGETRS ---
#define WRAP_DGETRS(F,A1,A2,A3,A4,A5,A6,A7,A8,INFO)   \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5,A6,A7,A8)

// --- DPOTRS ---
#define WRAP_DPOTRS(F,A1,A2,A3,A4,A5,A6,A7,INFO) \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5,A6,A7)


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
  fprintf(stderr, "Your lapack (atlas) version misses dgels function.\n");
#endif

// --- DGETRI ---
#define WRAP_DGETRI(F,A1,A2,A3,A4,INFO)         \
  INFO = F(CblasColMajor,A1,A2,A3,A4)

// --- DPOTRF ---
#define WRAP_DPOTRF(F,A1,A2,A3,A4,INFO)  \
  INFO = F(CblasColMajor,A1,A2,A3,A4)

// --- DGETRF ---
#define WRAP_DGETRF(F,A1,A2,A3,A4,A5,INFO)  \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5)


// --- DTRTRS ---
#if defined(HAS_LAPACK_dtrtrs)
#define WRAP_DTRTRS(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,INFO)  \
  INFO = F(CblasColMajor,A1,A2,A3,A4,A5,A6,A7,A8,A9)
#else
#define WRAP_DTRTRS(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,INFO)                \
  fprintf(stderr, "Your lapack(atlas) version misses dtrtrs function.\n");
#endif

#endif // SiconosAtlasLAPACK_H
