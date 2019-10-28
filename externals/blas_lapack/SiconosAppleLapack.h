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

#ifndef SiconosAppleLAPACK_H
#define SiconosAppleLAPACK_H

#include "SiconosBlas.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>


#define FCAST(T,X) (T *) (& X)
#define FCASTP(T,X) (T *) X

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
#define WRAP_DGESVD(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,INFO)      \
  __CLPK_integer lwork = -1;                                            \
  __CLPK_doublereal* work = (__CLPK_doublereal*)malloc(sizeof(*work));  \
  F(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,work,&lwork,INFO);               \
  lwork = (__CLPK_integer)work[0];                                      \
  work = realloc(work,lwork*sizeof(*work));                             \
  F(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,work,&lwork,INFO);               \
  free(work);                                                           \

// --- DGETRS ---
#define WRAP_DGETRS(F,A1,A2,A3,A4,A5,A6,A7,A8,INFO) \
  F(A1,A2,A3,A4,A5,A6,A7,A8,INFO)

// --- DPOTRS ---
#define WRAP_DPOTRS(F,A1,A2,A3,A4,A5,A6,A7,INFO) \
  F(A1,A2,A3,A4,A5,A6,A7,INFO)

// --- DGESV ---
#define WRAP_DGESV(F,A1,A2,A3,A4,A5,A6,A7,INFO) \
  F(A1,A2,A3,A4,A5,A6,A7,INFO)

// --- DPOSV ---
#define WRAP_DPOSV(F,A1,A2,A3,A4,A5,A6,A7,INFO) \
  F(A1,A2,A3,A4,A5,A6,A7,INFO)

// --- DGELS ---
#define WRAP_DGELS(F,A1,A2,A3,A4,A5,A6,A7,A8,INFO)                      \
  __CLPK_integer lwork = -1;                                            \
  __CLPK_doublereal* work = (__CLPK_doublereal*)malloc(sizeof(*work));  \
  F(A1,A2,A3,A4,A5,A6,A7,A8,work,&lwork,INFO);                          \
  lwork = (__CLPK_integer)work[0];                                      \
  work = realloc(work,lwork*sizeof(*work));                             \
  F(A1,A2,A3,A4,A5,A6,A7,A8,work,&lwork,INFO);                          \
  free(work);                                                           \

// --- DGETRI ---
#define WRAP_DGETRI(F,A1,A2,A3,A4,INFO)                                 \
  __CLPK_integer lwork = -1;                                            \
  __CLPK_doublereal* C_WORK;                                            \
  C_WORK = (__CLPK_doublereal*)malloc(sizeof *C_WORK);                  \
  assert(C_WORK);                                                       \
  F(A1,A2,A3,A4,C_WORK,&lwork,INFO);                                    \
  lwork = (__CLPK_integer) (C_WORK[0]);                                 \
  C_WORK = (__CLPK_doublereal*)realloc(C_WORK, lwork * sizeof *C_WORK); \
  F(A1,A2,A3,A4,C_WORK,&lwork,INFO);                                    \
  free(C_WORK);                                                         \

// --- DPOTRF ---
#define WRAP_DPOTRF(F,A1,A2,A3,A4,INFO)                                 \
  F(A1,A2,A3,A4,INFO);                                                  \

// --- DGETRF ---
#define WRAP_DGETRF(F,A1,A2,A3,A4,A5,INFO)      \
  F(A1,A2,A3,A4,A5,INFO)

// --- DTRTRS ---
#define WRAP_DTRTRS(F,A1,A2,A3,A4,A5,A6,A7,A8,A9,INFO) \
  F(A1,A2,A3,A4,A5,A6,A7,A8,A9,INFO)


#endif // SICONOSAPPLELAPACK_H
