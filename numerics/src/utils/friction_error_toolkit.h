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

#ifndef FRICTION_ERROR_TOOLKIT_H
#define FRICTION_ERROR_TOOLKIT_H

/*!\file friction_error_toolkit.h
 * \brief Shared helpers for friction contact error computations.
 */

#include "NumericsMatrix.h"  // for NumericsMatrix, NM_gemv, NM_tgemv
#include "SiconosBlas.h"     // for cblas_dcopy, cblas_dnrm2, cblas_daxpy
#include "sanitizer.h"       // for cblas_dcopy_msan

/** Recompute the local velocity from the global one.
 *
 * Computes velocity = b + H^T * globalVelocity.
 *
 * \param[in] H localization matrix
 * \param[in] b local offset vector (size m)
 * \param[in] globalVelocity global velocity vector (size n)
 * \param[in] m dimension of the local velocity vector
 * \param[in,out] velocity local velocity vector (size m), overwritten on output
 */
static inline void friction_error_toolkit_recompute_local_velocity(
    NumericsMatrix* H, const double* b, const double* globalVelocity, int m,
    double* velocity) {
  cblas_dcopy(m, b, 1, velocity, 1);
  NM_tgemv(1, H, globalVelocity, 1.0, velocity);
}

/** Compute the global primal residual for a global friction contact problem.
 *
 * Computes rhs = -M * globalVelocity + H * reaction + q and returns its Euclidean norm.
 * Optionally returns the norms of H * reaction and M * globalVelocity, which are needed
 * for relative scaling.  If has_reaction is zero the H * reaction term is skipped
 * (used for problems with no contacts).
 *
 * \param[in] M global mass/dynamics matrix
 * \param[in] H localization matrix
 * \param[in] q global right-hand side vector (size n)
 * \param[in] reaction local reaction vector (size m)
 * \param[in] globalVelocity global velocity vector (size n)
 * \param[in] n global dimension
 * \param[in] has_reaction if non-zero, include H * reaction; otherwise skip it
 * \param[in,out] work workspace vector of size n (overwritten)
 * \param[in,out] rhs output residual vector of size n (overwritten)
 * \param[out] norm_Hr optional output for ||H * reaction|| (may be NULL)
 * \param[out] norm_Mv optional output for ||M * globalVelocity|| (may be NULL)
 * \return Euclidean norm of the residual vector rhs
 */
static inline double friction_error_toolkit_global_primal_residual(
    NumericsMatrix* M, NumericsMatrix* H, const double* q, const double* reaction,
    const double* globalVelocity, size_t n, int has_reaction, double* work, double* rhs,
    double* norm_Hr, double* norm_Mv) {
  const int ni = (int)n;
  cblas_dcopy_msan(ni, q, 1, rhs, 1);

  if (has_reaction) {
    NM_gemv(1.0, H, reaction, 0.0, work);
    const double nHr = cblas_dnrm2(ni, work, 1);
    cblas_daxpy(ni, 1.0, work, 1, rhs, 1);
    if (norm_Hr) *norm_Hr = nHr;
  } else {
    if (norm_Hr) *norm_Hr = 0.0;
  }

  NM_gemv(-1.0, M, globalVelocity, 0.0, work);
  const double nMv = cblas_dnrm2(ni, work, 1);
  cblas_daxpy(ni, 1.0, work, 1, rhs, 1);

  if (norm_Mv) *norm_Mv = nMv;

  return cblas_dnrm2(ni, rhs, 1);
}

/** Compute the global dual residual for a global friction contact problem.
 *
 * Computes rhs = H^T * globalVelocity + b - velocity and returns its Euclidean norm.
 * Optionally returns the norm of H^T * globalVelocity for relative scaling.
 *
 * \param[in] H localization matrix
 * \param[in] b local offset vector (size m)
 * \param[in] globalVelocity global velocity vector (size n)
 * \param[in] velocity local velocity vector (size m)
 * \param[in] m local dimension
 * \param[in,out] work workspace vector of size m (overwritten)
 * \param[in,out] rhs output residual vector of size m (overwritten)
 * \param[out] norm_Htv optional output for ||H^T * globalVelocity|| (may be NULL)
 * \return Euclidean norm of the residual vector rhs
 */
static inline double friction_error_toolkit_global_dual_residual(
    NumericsMatrix* H, const double* b, const double* globalVelocity, const double* velocity,
    size_t m, double* work, double* rhs, double* norm_Htv) {
  const int mi = (int)m;
  cblas_dcopy_msan(mi, b, 1, rhs, 1);

  NM_tgemv(1, H, globalVelocity, 0.0, work);
  const double nHtv = cblas_dnrm2(mi, work, 1);
  cblas_daxpy(mi, 1.0, work, 1, rhs, 1);

  cblas_daxpy(mi, -1.0, velocity, 1, rhs, 1);

  if (norm_Htv) *norm_Htv = nHtv;

  return cblas_dnrm2(mi, rhs, 1);
}

#endif /* FRICTION_ERROR_TOOLKIT_H */
