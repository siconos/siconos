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


