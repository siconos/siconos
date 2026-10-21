/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
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
#include "fc3d_compute_rho.h"  // for fc3d_compute_rho

#include <assert.h>  // for assert
#include <math.h>    // for fabs
#include <stdlib.h>  // for calloc, free

#include "FrictionContactProblem.h"  // for frictionContactProblem_new_with_data, ...
#include "NumericsMatrix.h"          // for NM_create_with_data, NM_DENSE

int main() {
  int nc = 2;
  int dim = 3;
  int n = dim * nc;

  double* Mdata = (double*)calloc((size_t)n * n, sizeof(double));
  for (int i = 0; i < n; ++i) Mdata[i + i * n] = 1.0;
  Mdata[1 + 1 * n] = 4.0;
  Mdata[2 + 2 * n] = 9.0;
  Mdata[4 + 4 * n] = 2.0;
  Mdata[5 + 5 * n] = 5.0;

  double* q = (double*)calloc((size_t)n, sizeof(double));
  double* mu = (double*)calloc((size_t)nc, sizeof(double));
  for (int i = 0; i < nc; ++i) mu[i] = 1.0;

  NumericsMatrix* M = NM_create_with_data(NM_DENSE, n, n, Mdata);
  FrictionContactProblem* problem =
      frictionContactProblem_new_with_data(dim, (size_t)nc, M, q, mu);

  double* rho = (double*)calloc((size_t)3 * nc, sizeof(double));

  /* Constant strategy */
  int info = fc3d_compute_rho(
      problem, SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT, 0.5, rho);
  assert(info == 0);
  for (int i = 0; i < 3 * nc; ++i) assert(fabs(rho[i] - 0.5) < 1e-12);

  /* Split spectral norm strategy */
  info = fc3d_compute_rho(
      problem, SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM, 0.0, rho);
  assert(info == 0);
  assert(fabs(rho[0] - 1.0) < 1e-12);
  assert(rho[1] > 0.0);
  assert(rho[2] > 0.0);
  assert(fabs(rho[3] - 1.0) < 1e-12);
  assert(rho[4] > 0.0);
  assert(rho[5] > 0.0);

  /* Split spectral norm cond strategy */
  info = fc3d_compute_rho(
      problem, SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND, 0.0,
      rho);
  assert(info == 0);
  assert(fabs(rho[0] - 1.0) < 1e-12);
  assert(rho[1] > 0.0);
  assert(rho[2] > 0.0);
  assert(fabs(rho[3] - 1.0) < 1e-12);
  assert(rho[4] > 0.0);
  assert(rho[5] > 0.0);

  /* Spectral norm strategy */
  info = fc3d_compute_rho(
      problem, SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM, 0.0, rho);
  assert(info == 0);
  assert(fabs(rho[0] - 1.0) < 1e-12);
  assert(rho[1] > 0.0);
  assert(rho[2] > 0.0);

  free(rho);
  frictionContactProblem_free(problem);

  return 0;
}
