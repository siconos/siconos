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

#include "gfc2d_VI_ExtraGradient.h"

#include <stdio.h>   // for printf
#include <stdlib.h>  // for free, malloc, calloc
#include <string.h>  // for memcpy

#include "GlobalFrictionContactProblem.h"        // for GlobalFrictionContac...
#include "GlobalFrictionContactProblem_as_VI.h"  // for GlobalFrictionContac...
#include "NumericsFwd.h"                         // for VariationalInequality
#include "NumericsMatrix.h"                      // for NumericsMatrix
#include "SiconosBlas.h"                         // for cblas_dnrm2
#include "SolverOptions.h"                       // for SolverOptions, SICON...
#include "VariationalInequality.h"               // for VariationalInequality
#include "VariationalInequality_Solvers.h"       // for variationalInequalit...
#include "fc3d_short_names.h"                    // for GFC2D_VI_EG
#include "gfc2d_compute_error.h"                 // for gfc2d_compute_error
#include "numerics_verbose.h"
#include "siconos_debug.h"  // for DEBUG_EXPR, DEBUG_BEGIN

/* Solver registration system */
#include "numerics_errors.h"
#include "solver_registry.h"

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

void gfc2d_VI_ExtraGradient(GlobalFrictionContactProblem* problem, double* reaction,
                            double* velocity, double* globalVelocity, int* info,
                            SolverOptions* options) {
  DEBUG_BEGIN("gfc2d_VI_ExtraGradient(GlobalFrictionContactProblem* problem, ... \n");
  DEBUG_EXPR(verbose = 1;);
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  /* Dimension of the local problem */
  int m = 2 * nc;
  int n = problem->M->size0;

  DEBUG_EXPR(NM_vector_display(reaction, m););
  DEBUG_EXPR(NM_vector_display(velocity, m););
  DEBUG_EXPR(NM_vector_display(globalVelocity, n););

  VariationalInequality* vi = (VariationalInequality*)malloc(sizeof(VariationalInequality));

  vi->F = &Function_VI_GFC2D;
  vi->ProjectionOnX = &Projection_VI_GFC2D;

  double error = 1e24;

  GlobalFrictionContactProblem_as_VI* gfc2d_as_vi =
      (GlobalFrictionContactProblem_as_VI*)malloc(sizeof(GlobalFrictionContactProblem_as_VI));
  vi->env = gfc2d_as_vi;
  vi->size = n + m;

  /* Set the norm of the VI to the norm of problem->q */
  vi->normVI = cblas_dnrm2(n, problem->q, 1);
  vi->istheNormVIset = 1;

  gfc2d_as_vi->vi = vi;
  gfc2d_as_vi->gfc3d = problem;

  DEBUG_EXPR(NV_display(reaction, m););
  DEBUG_EXPR(NV_display(globalVelocity, n););
  double* z = (double*)malloc((n + m) * sizeof(double));
  double* Fz = (double*)calloc((n + m), sizeof(double));

  memcpy(z, globalVelocity, n * sizeof(double));
  memcpy(&z[n], reaction, m * sizeof(double));
  DEBUG_EXPR(NV_display(z, m + n););
  DEBUG_EXPR(NV_display(Fz, m + n););

  variationalInequality_ExtraGradient(vi, z, Fz, info, options);

  memcpy(globalVelocity, z, n * sizeof(double));
  memcpy(reaction, &z[n], m * sizeof(double));
  memcpy(velocity, &Fz[n], m * sizeof(double));

  free(z);
  free(Fz);

  double norm_q = cblas_dnrm2(n, problem->q, 1);
  double norm_b = cblas_dnrm2(m, problem->b, 1);

  /* **** Criterium convergence **** */
  gfc2d_compute_error(problem, reaction, velocity, globalVelocity,
                      options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);

  DEBUG_EXPR(NM_vector_display(reaction, m));
  DEBUG_EXPR(NM_vector_display(velocity, m));
  DEBUG_EXPR(NM_vector_display(globalVelocity, n));

  if (verbose > 0) {
    printf(
        "--------------- GFC2D - VI Extra Gradient (VI_EG) - #Iteration %i Final Residual = "
        "%14.7e\n",
        options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  }

  free(vi);
  free(gfc2d_as_vi);

  DEBUG_END("gfc2d_VI_ExtraGradient(GlobalFrictionContactProblem* problem, ... \n")
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers GFC2D_VI_EG in the global solver registry, enabling:
 * - Dynamic solver lookup by ID
 * - Runtime solver introspection
 * - Elimination of giant switch statements in drivers
 */

void gfc2d_VI_extra_gradient_set_default(SolverOptions* options) {
  /* VI_EG doesn't use internal solvers, but we call the VI set_default for proper
   * initialization */
  variationalInequality_ExtraGradient_set_default(options);
}

static int gfc2d_vi_eg_init_wrap(void* problem, SolverOptions* options) {
  /* set_default already called by solver_options_create */
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int gfc2d_vi_eg_solve_wrap(void* problem, double* reaction, double* velocity,
                                  double* globalVelocity, SolverOptions* options) {
  int info = NUMERICS_OK;
  gfc2d_VI_ExtraGradient((GlobalFrictionContactProblem*)problem, reaction, velocity,
                         globalVelocity, &info, options);
  return info;
}

static void gfc2d_vi_eg_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER_3VAR(GFC2D_VI_EG, "GFC2D_VI_EG",
                     "Variational Inequality Extra Gradient for 2D Global Friction Contact",
                     gfc2d_vi_eg_init_wrap, gfc2d_vi_eg_solve_wrap, gfc2d_vi_eg_free_wrap,
                     NULL,                                /* error function */
                     gfc2d_VI_extra_gradient_set_default, /* set_default */
                     1000,                                /* default_max_iter */
                     1e-4,                                /* default_tol */
                     0 /* is_local_solver */);
