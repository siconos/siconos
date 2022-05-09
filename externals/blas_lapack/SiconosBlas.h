/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#ifndef SiconosBlas_H
#define SiconosBlas_H
#include "SiconosConfig.h"

#if defined(__cplusplus)
extern "C"
{
#endif

// tells include-what-you-use to keep this file
// and not to suggest cblas.h or alike.
// IWYU pragma: begin_exports
#if defined(HAS_MKL_CBLAS)
#include <mkl_cblas.h>
#elif defined(HAS_MATLAB_BLAS)
#include <blas.h>
#define cblas_daxpy daxpy
#define cblas_dcopy dcopy
#define cblas_ddot ddot
#define cblas_dgemm dgemm
#define cblas_dgemv dgemv
#define cblas_dnrm2 dnrm2
#define cblas_dscal dscal
#else
#include <cblas.h>
#endif
// IWYU pragma: end_exports

#ifdef __cplusplus
}
#undef restrict
#define restrict __restrict
#endif


static inline double* NMD_row_rmajor(double* restrict mat, unsigned ncols, unsigned rindx)
{
  return &mat[rindx*ncols];
}

static inline void NMD_copycol_rmajor(int nrows, double* col, double* restrict mat, int ncols, unsigned cindx)
{
  cblas_dcopy(nrows, col, 1, &mat[cindx], ncols);
}

static inline void NMD_dense_gemv(int nrows, int ncols, double alpha, double* restrict mat, double* restrict y, double beta, double* restrict x)
{
  cblas_dgemv(CblasColMajor, CblasTrans, ncols, nrows, alpha, mat, ncols, y, 1, beta, x, 1);
}

#endif // SiconosBlas_H


