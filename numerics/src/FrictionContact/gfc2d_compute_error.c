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

#include "gfc2d_compute_error.h"

#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for fabs, sqrt
#include <stdlib.h>  // for NULL, realloc

#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "fc2d_compute_error.h" 
#include "NumericsMatrix.h"                // for NM_gemv
#include "SiconosBlas.h"                   // for cblas_dcopy, cblas_dnrm2, cblas_daxpy
#include "SolverOptions.h"                 // for SolverOptions
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "sanitizer.h"  // for cblas_dcopy_msan

//#define DEBUG_MESSAGES
#include "siconos_debug.h"
#include "numerics_errors.h"

#define MIN_RELATIVE_SCALING sqrt(DBL_EPSILON)

int gfc2d_compute_error(GlobalFrictionContactProblem* problem, double* reaction,
                        double* velocity, double* globalVelocity, double tolerance,
                        SolverOptions* options, double norm_q, double norm_b,
                        double* error) {
  (void)velocity;
  (void)tolerance;
  (void)norm_b;

  CHECK_NULL(problem);
  CHECK_NULL(reaction);
  CHECK_NULL(globalVelocity);
  CHECK_NULL(error);
  CHECK_MATRIX(problem->M);
  CHECK_MATRIX(problem->H);
  CHECK_NULL(problem->q);

  const size_t n = problem->M->size0;
  const int nc = problem->numberOfContacts;
  const int dim = problem->dimension;
  const int m = nc * dim;
  if (!options->dWork || options->dWorkSize < 2 * n) {
    double* p = (double*)realloc(options->dWork, 2 * n * sizeof(double));
    if (!p) {
      numerics_error_log("gfc2d_compute_error", "bad alloc");
      return NUMERICS_ERR_MEMORY;
    }
    options->dWork = p;
    options->dWorkSize = 2 * n;
  }

  double* tmp = options->dWork;
  double* rhs = &options->dWork[n];

  /* rhs = q */
  cblas_dcopy_msan(n, problem->q, 1, rhs, 1);

  /* tmp = H * reaction; rhs += tmp */
  if (problem->numberOfContacts > 0) {
    NM_gemv(1.0, problem->H, reaction, 0.0, tmp);
    cblas_daxpy(n, 1.0, tmp, 1, rhs, 1);
  }
  const double norm_Hr = cblas_dnrm2(n, tmp, 1);

  /* tmp = -M * globalVelocity; rhs += tmp */
  NM_gemv(-1.0, problem->M, globalVelocity, 0.0, tmp);
  cblas_daxpy(n, 1.0, tmp, 1, rhs, 1);

  const double norm_Mv = cblas_dnrm2(n, tmp, 1);

  const double relative_scaling = fmax(norm_q, fmax(norm_Mv, norm_Hr));
  *error = cblas_dnrm2(n, rhs, 1);
  DEBUG_PRINTF("absolute error  of -M v + H R + q = %e\n", *error);
  
 /* we re-compute local velocity */
  /* the error in the equation u = H^T v +b is then accurate at the machine precision */

  cblas_dcopy(m, problem->b, 1, velocity, 1);
  NM_tgemv(1, problem->H, globalVelocity, 1.0, velocity);
  double norm_u = cblas_dnrm2(m, velocity, 1);
  DEBUG_PRINTF("norm of velocity %e\n", norm_u);


  double error_complementarity = 0.0;
  double worktmp[2];
  for (int ic = 0; ic < nc; ic++) {
    fc2d_unitary_compute_and_add_error(&reaction[ic * dim], &velocity[ic * dim],
                                       problem->mu[ic],
                                       &error_complementarity, worktmp);
  }

  error_complementarity = sqrt(error_complementarity);

  DEBUG_PRINTF("absolute error in cone  complementarity= %e\n", error_complementarity);

  *error += error_complementarity;
  DEBUG_PRINTF("gfc2d_compute_error: absolute error = %e\n", *error);
  if (fabs(relative_scaling) > MIN_RELATIVE_SCALING) {
    *error = *error / relative_scaling;
  }

  DEBUG_PRINTF("gfc2d_compute_error: relative error = %e\n", *error);

  return NUMERICS_OK;
}
