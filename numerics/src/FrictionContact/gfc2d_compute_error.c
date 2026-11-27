// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "gfc2d_compute_error.h"

#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for fabs, sqrt
#include <stdlib.h>  // for NULL, realloc

#include "FrictionContact_options.h"
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NumericsMatrix.h"                // for NM_gemv
#include "SiconosBlas.h"                   // for cblas_dcopy, cblas_dnrm2, cblas_daxpy
#include "SolverOptions.h"                 // for SolverOptions
#include "error_computation.h"             // for error_compute_alloc_pair
#include "fc2d_compute_error.h"
#include "friction_error_toolkit.h"  // for friction_error_toolkit_...
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "sanitizer.h"  // for cblas_dcopy_msan

// #define DEBUG_MESSAGES
#include "numerics_errors.h"
#include "siconos_debug.h"

#define MIN_RELATIVE_SCALING sqrt(DBL_EPSILON)

int gfc2d_compute_error(GlobalFrictionContactProblem* problem, double* reaction,
                        double* velocity, double* globalVelocity, double tolerance,
                        SolverOptions* options, double norm_q, double norm_b, double* error) {
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
  double* mu = problem->mu;

  const NaturalMapErrorScaling defaults = {SICONOS_FRICTION_2D_SCALING_EUCLIDEAN,
                                           SICONOS_FRICTION_2D_SCALING_SCOPE_PER_CONE,
                                           SICONOS_FRICTION_2D_SCALING_REDUCTION_MAX, 1.0};
  NaturalMapErrorScaling sc = natural_map_error_scaling(options, defaults);

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

  double norm_Hr = 0.0, norm_Mv = 0.0;
  *error = friction_error_toolkit_global_primal_residual(problem->M, problem->H, problem->q,
                                                         reaction, globalVelocity, n, nc > 0,
                                                         tmp, rhs, &norm_Hr, &norm_Mv);

  const double relative_scaling = fmax(norm_q, fmax(norm_Mv, norm_Hr));
  DEBUG_PRINTF("absolute error  of -M v + H R + q = %e\n", *error);

  /* we re-compute local velocity */
  /* the error in the equation u = H^T v +b is then accurate at the machine precision */

  friction_error_toolkit_recompute_local_velocity(problem->H, problem->b, globalVelocity, m,
                                                  velocity);
  DEBUG_PRINTF("norm of velocity %e\n", cblas_dnrm2(m, velocity, 1));

  double error_complementarity = 0.0;
  if (nc > 0) {
    /* natural-map error on the 2D friction cones (absolute if the norm is NO) */
    double* comp_work = NULL;
    double* comp_tmp = NULL;
    if (error_compute_alloc_pair(nc, m, &comp_work, &comp_tmp)) {
      return numerics_error("gfc2d_compute_error", "bad alloc");
    }
    error_complementarity =
        fc2d_compute_complementarity_error(nc, reaction, velocity, mu, sc.norm, sc.scope,
                                           sc.reduction, sc.alpha, comp_work, comp_tmp);
    free(comp_work);
    free(comp_tmp);
  }

  DEBUG_PRINTF("absolute error in cone  complementarity= %e\n", error_complementarity);

  *error += error_complementarity;
  DEBUG_PRINTF("gfc2d_compute_error: absolute error = %e\n", *error);
  if (fabs(relative_scaling) > MIN_RELATIVE_SCALING) {
    *error = *error / relative_scaling;
  }

  DEBUG_PRINTF("gfc2d_compute_error: relative error = %e\n", *error);

  return NUMERICS_OK;
}
