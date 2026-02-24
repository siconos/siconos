
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
#include <assert.h>  // for assert
#include <math.h>    // for sqrt
#include <stdio.h>   // for fprintf, stderr
#include <stdlib.h>  // for malloc, free

#include "ConvexQP.h"                                  // for ConvexQP
#include "ConvexQP_Solvers.h"                          // for convexQP_ADMM_...
#include "ConvexQP_cst.h"                              // for SICONOS_CONVEX...
#include "FrictionContact_options.h"                              // for SICONOS_FRICTI...
#include "GlobalFrictionContactProblem.h"              // for GlobalFriction...
#include "GlobalFrictionContactProblem_as_ConvexQP.h"  // for GlobalFriction...
#include "NumericsFwd.h"                               // for ConvexQP, Solv...
#include "NumericsMatrix.h"                            // for NM_clear, NM_tr...
#include "SiconosBlas.h"                               // for cblas_dcopy
#include "SolverOptions.h"                             // for SolverOptions
#include "gfc3d_Solvers.h"                             // for gfc3d_set_inte...
#include "gfc3d_compute_error.h"                       // for gfc3d_compute_...
#include "fc3d_short_names.h"                          // for GFC3D_ACLMFP
#include "numerics_verbose.h"
#include "siconos_debug.h"                             // for DEBUG_EXPR

/* Solver registration system */
#include "utils/solver_registry.h"
#include "utils/numerics_errors.h"

/** pointer to function used to call internal solver for proximal point solver */
typedef void (*internalSolverPtr)(ConvexQP *, double *, double *, double *, double *, int *,
                                  SolverOptions *);
void gfc3d_ACLMFixedPoint(GlobalFrictionContactProblem *restrict problem,
                          double *restrict reaction, double *restrict velocity,
                          double *restrict globalVelocity, int *restrict info,
                          SolverOptions *restrict options) {
  /* verbose=1; */

  /* int and double parameters */
  int *iparam = options->iparam;
  double *dparam = options->dparam;

  /* Number of contacts */
  int nc = problem->numberOfContacts;
  int n = problem->M->size0;
  int m = 3 * nc;

  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  double norm_q = cblas_dnrm2(n, problem->q, 1);
  double norm_b = cblas_dnrm2(m, problem->b, 1);

  if (options->numberOfInternalSolvers < 1) {
    numerics_error("gfc3d_ACLMFixedpoint",
                   "The ACLM Fixed Point method needs options for the internal solvers, "
                   "options[0].numberOfInternalSolvers should be >1");
  }

  SolverOptions *internalsolver_options = options->internalSolvers[0];

  if (verbose > 0) {
    solver_options_print(options);
  }

  /*****  Fixed Point Iterations *****/
  int iter = 0;      /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  internalSolverPtr internalsolver;

  ConvexQP *cqp = (ConvexQP *)malloc(sizeof(ConvexQP));
  cqp->size = n;
  cqp->m = m;
  cqp->M = problem->M;

  cqp->q = (double *)malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) {
    cqp->q[i] = -problem->q[i];
  }
  DEBUG_EXPR(NM_display(problem->H));

  cqp->A = NM_transpose(problem->H);
  DEBUG_EXPR(NM_display(cqp->A));
  cqp->b = (double *)malloc(m * sizeof(double));
  cqp->ProjectionOnC = &Projection_ConvexQP_GFC3D_DualCone;
  cqp->normConvexQP = norm_q;
  cqp->istheNormConvexQPset = 1;
  double *w = (double *)malloc(n * sizeof(double));

  GlobalFrictionContactProblem_as_ConvexQP *gfc3d_as_cqp =
      (GlobalFrictionContactProblem_as_ConvexQP *)malloc(
          sizeof(GlobalFrictionContactProblem_as_ConvexQP));
  cqp->env = gfc3d_as_cqp;
  gfc3d_as_cqp->cqp = cqp;
  gfc3d_as_cqp->gfc3d = problem;
  gfc3d_as_cqp->options = options;
  if (internalsolver_options->solverId == SICONOS_CONVEXQP_ADMM) {
    numerics_printf_verbose(1,
                            " ========================== set ADMM solver internal ConveQP "
                            "problem ==========================\n");
    internalsolver = &convexQP_ADMM;
    convexQP_ADMM_init(cqp, options->internalSolvers[0]);
  } else {
    fprintf(stderr, "Numerics, gfc3d_ACLMFixedPoint failed. Unknown internal solver.\n");
    exit(EXIT_FAILURE);
  }

  double normUT;
  int cumul_iter = 0;
  while ((iter < itermax) && (hasNotConverged > 0)) {
    ++iter;
    // internal solver for the regularized problem

    /* Compute the value of the initial value of b */
    cblas_dcopy(m, problem->b, 1, cqp->b, 1);
    for (int ic = 0; ic < nc; ic++) {
      normUT = sqrt(velocity[ic * 3 + 1] * velocity[ic * 3 + 1] +
                    velocity[ic * 3 + 2] * velocity[ic * 3 + 2]);
      cqp->b[3 * ic] += problem->mu[ic] * normUT;
    }

    gfc3d_set_internalsolver_tolerance(problem, options, internalsolver_options, error);

    (*internalsolver)(cqp, globalVelocity, w, reaction, velocity, info,
                      internalsolver_options);

    cumul_iter += internalsolver_options->iparam[SICONOS_IPARAM_ITER_DONE];
    /* **** Criterium convergence **** */

    gfc3d_compute_error(problem, reaction, velocity, globalVelocity, tolerance, options,
                        norm_q, norm_b, &error);

    numerics_printf_verbose(1, "---- GFC3D - ACLMFP - Iteration %i Residual = %14.7e", iter,
                            error);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }

  numerics_printf_verbose(1, "---- GFC3D - ACLMFP - # Iteration %i Final Residual = %14.7e",
                          iter, error);
  numerics_printf_verbose(1, "---- GFC3D - ACLMFP - #              internal iteration = %i",
                          cumul_iter);

  NM_clear(cqp->A);
  free(cqp->b);
  free(cqp->q);
  free(w);

  if (internalsolver_options->solverId == SICONOS_CONVEXQP_ADMM) {
    convexQP_ADMM_free(cqp, options->internalSolvers[0]);
  }
  free(cqp);

  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;
}

void gfc3d_aclmfp_set_default(SolverOptions *options) {
  options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] =
      SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE;
  options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] = 2.0;

  // Internal solver - allocate if needed
  if (options->numberOfInternalSolvers == 0) {
    options->numberOfInternalSolvers = 1;
    options->internalSolvers = calloc(1, sizeof(SolverOptions*));
  }
  assert(options->numberOfInternalSolvers == 1);
  options->internalSolvers[0] = solver_options_create(SICONOS_CONVEXQP_ADMM);
  options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers GFC3D_ACLMFP in the global solver registry, enabling:
 * - Dynamic solver lookup by ID
 * - Runtime solver introspection
 * - Elimination of giant switch statements in drivers
 */

static int gfc3d_aclmfp_init_wrap(void* problem, SolverOptions* options) {
  gfc3d_aclmfp_set_default(options);
  return NUMERICS_OK;
}

static int gfc3d_aclmfp_solve_wrap(void* problem, double* reaction,
                                   double* velocity, SolverOptions* options) {
  int info = NUMERICS_OK;
  // For global solvers, we need to handle globalVelocity separately
  GlobalFrictionContactProblem* gfc3d_problem = (GlobalFrictionContactProblem*)problem;
  int n = gfc3d_problem->M->size0;
  double* globalVelocity = (double*)calloc(n, sizeof(double));
  gfc3d_ACLMFixedPoint(gfc3d_problem, reaction, velocity, globalVelocity, &info, options);
  free(globalVelocity);
  return info;
}

static void gfc3d_aclmfp_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER_WITH_DEFAULT(GFC3D_ACLMFP, "GFC3D_ACLMFP",
                "Alart-Curnier Lemke Fixed Point for 3D Global Friction Contact",
                gfc3d_aclmfp_init_wrap,
                gfc3d_aclmfp_solve_wrap,
                gfc3d_aclmfp_free_wrap,
                NULL,  /* error function */
                gfc3d_aclmfp_set_default,  /* set_default */
                1000,  /* default_max_iter */
                1e-4,  /* default_tol */
                0      /* is_local_solver */);
