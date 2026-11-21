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
#include "nonsmooth_newton_levenberg_marquardt.h"

#include <float.h>  // for DBL_MIN

#include "NumericsMatrix.h"  // for NM_transpose, NM_multiply, NM_add, NM_eye, NM_LU_solve

int nonsmooth_newton_lm_is_stalled(double residual_norm, double previous_residual_norm,
                                   int* stalled_iterations) {
  if (residual_norm > NONSMOOTH_NEWTON_LM_STALL_RATIO * previous_residual_norm)
    ++(*stalled_iterations);
  else
    *stalled_iterations = 0;
  return *stalled_iterations >= NONSMOOTH_NEWTON_LM_STALL_ITERATIONS;
}

int nonsmooth_newton_lm_direction(NumericsMatrix* jacobian, const double* residual,
                                  double lambda, double* direction) {
  /* The normal matrix is built with sparse products: work on a sparse copy of J when it is
   * stored as a dense or a sparse-block matrix (the exact zeros are dropped). */
  NumericsMatrix* J = jacobian;
  if (jacobian->storageType != NM_SPARSE) J = NM_copy_to_sparse(jacobian, DBL_MIN);

  NumericsMatrix* jacobian_t = NM_transpose(J);
  NumericsMatrix* jtj = NM_multiply(jacobian_t, J);
  NumericsMatrix* eye = NM_eye(J->size1);
  NumericsMatrix* lm_matrix = NM_add(lambda, eye, 1.0, jtj);

  NM_gemv(-1.0, jacobian_t, residual, 0.0, direction);  // -J^T Psi
  const int info = NM_LU_solve(NM_preserve(lm_matrix), direction, 1);

  NM_free(lm_matrix);
  NM_free(eye);
  NM_free(jtj);
  NM_free(jacobian_t);
  if (J != jacobian) NM_free(J);
  return info;
}
