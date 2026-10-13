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

#include "gfc2d_VI_FixedPointProjection.h"

#include <math.h>    // for fabs
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
#include "fc3d_short_names.h"                    // for GFC2D_VI_FPP
#include "gfc2d_compute_error.h"                 // for gfc2d_compute_error
#include "numerics_verbose.h"
#include "siconos_debug.h"  // for DEBUG_EXPR, DEBUG_BEGIN

/* Solver registration system */
#include "numerics_errors.h"
#include "solver_registry.h"

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

/**
 * \brief Classify a single 2D contact into a discrete status.
 *
 * \return 0 = no-contact, 1 = stick, 2 = slip in the positive tangent
 *         direction, 3 = slip in the negative tangent direction.
 */
static int gfc2d_contact_status_2d(double rN, double rT, double uT, double mu) {
  const double tol = 1e-10;
  const double radius = mu * rN;
  if (rN <= tol) {
    return 0; /* no-contact */
  }
  if (fabs(rT) < radius - tol && fabs(uT) <= tol) {
    return 1; /* stick */
  }
  if (uT > tol) {
    return 2; /* slip+ */
  }
  if (uT < -tol) {
    return 3; /* slip- */
  }
  /* Near the boundary: use the sign of the friction force. */
  return (rT < 0.0) ? 2 : 3;
}

typedef struct {
  GlobalFrictionContactProblem* problem;
  int n_global;
  int iteration;
  int* previous_status;
  int* current_status;
} GFC2D_VIFPP_StatusEnv;

/**
 * \brief Per-iteration callback that prints the 2D contact status counts and
 *        the percentage of contacts whose status did not change since the
 *        previous VI-FPP iteration.
 */
static void gfc2d_vi_fpp_status_callback(void* env, size_t size, double* x, double* velocity,
                                         double error, void* extra_data) {
  (void)size;
  (void)error;
  (void)extra_data;

  GFC2D_VIFPP_StatusEnv* status_env = (GFC2D_VIFPP_StatusEnv*)env;
  GlobalFrictionContactProblem* problem = status_env->problem;
  const int nc = problem->numberOfContacts;
  const int n_global = status_env->n_global;
  const double* reaction = &x[n_global];
  const double* local_velocity = &velocity[n_global];

  unsigned int no_contact = 0, stick = 0, slip_right = 0, slip_left = 0;
  unsigned int stabilized = 0;

  for (int i = 0; i < nc; ++i) {
    const int status = gfc2d_contact_status_2d(reaction[2 * i], reaction[2 * i + 1],
                                               local_velocity[2 * i + 1], problem->mu[i]);
    status_env->current_status[i] = status;
    if (status == 0)
      no_contact++;
    else if (status == 1)
      stick++;
    else if (status == 2)
      slip_right++;
    else
      slip_left++;
    if (status_env->previous_status[i] != -1 && status_env->previous_status[i] == status) {
      stabilized++;
    }
    status_env->previous_status[i] = status;
  }

  const double stabilized_pct = (nc > 0) ? (100.0 * stabilized) / nc : 0.0;

  if (verbose > 0) {
    if (status_env->iteration == 1) {
      printf("Contact status: no-c | stick | slip+ | slip- | stab%%\n");
    }
    printf("                %4u | %5u | %5u | %5u | %5.1f%%\n", no_contact, stick, slip_right,
           slip_left, stabilized_pct);
  }
  status_env->iteration++;
}

void gfc2d_VI_FixedPointProjection(GlobalFrictionContactProblem* problem, double* reaction,
                                   double* velocity, double* globalVelocity, int* info,
                                   SolverOptions* options) {
  DEBUG_BEGIN("gfc2d_VI_FixedPointProjection(GlobalFrictionContactProblem* problem, ... \n");
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

  /* Set up a per-iteration callback that prints 2D contact statuses.
   * We preserve any user-defined callback that may already be installed. */
  GFC2D_VIFPP_StatusEnv status_env;
  status_env.problem = problem;
  status_env.n_global = n;
  status_env.iteration = 1;
  status_env.previous_status = (int*)malloc(nc * sizeof(int));
  status_env.current_status = (int*)malloc(nc * sizeof(int));
  for (int i = 0; i < nc; ++i) {
    status_env.previous_status[i] = -1;
  }

  Callback status_callback;
  status_callback.env = &status_env;
  status_callback.collectStatsIteration = gfc2d_vi_fpp_status_callback;

  Callback* saved_callback = options->callback;
  options->callback = &status_callback;

  variationalInequality_FixedPointProjection(vi, z, Fz, info, options);

  options->callback = saved_callback;
  free(status_env.previous_status);
  free(status_env.current_status);

  memcpy(globalVelocity, z, n * sizeof(double));
  memcpy(reaction, &z[n], m * sizeof(double));
  memcpy(velocity, &Fz[n], m * sizeof(double));

  free(z);
  free(Fz);

  /* **** Criterium convergence **** */
  double norm_q = cblas_dnrm2(n, problem->q, 1);
  double norm_b = cblas_dnrm2(m, problem->b, 1);

  gfc2d_compute_error(problem, reaction, velocity, globalVelocity,
                      options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);

  DEBUG_EXPR(NM_vector_display(reaction, m));
  DEBUG_EXPR(NM_vector_display(velocity, m));
  DEBUG_EXPR(NM_vector_display(globalVelocity, n));

  if (verbose > 0) {
    unsigned int no_contact = 0, stick = 0, slip_right = 0, slip_left = 0;
    for (int i = 0; i < nc; ++i) {
      const int status = gfc2d_contact_status_2d(reaction[2 * i], reaction[2 * i + 1],
                                                 velocity[2 * i + 1], problem->mu[i]);
      if (status == 0)
        no_contact++;
      else if (status == 1)
        stick++;
      else if (status == 2)
        slip_right++;
      else
        slip_left++;
    }
    printf(
        "---- GFC2D - VI_FPP - contact status: "
        "no-contact=%u, stick=%u, slip+=%u, slip-=%u\n",
        no_contact, stick, slip_right, slip_left);
    printf(
        "--------------- GFC2D - VI Fixed Point Projection (VI_FPP) - #Iteration %i Final "
        "Residual = %14.7e\n",
        options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  }
  free(vi);
  free(gfc2d_as_vi);

  DEBUG_END("gfc2d_VI_FixedPointProjection(GlobalFrictionContactProblem* problem, ... \n")
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers GFC2D_VI_FPP in the global solver registry, enabling:
 * - Dynamic solver lookup by ID
 * - Runtime solver introspection
 * - Elimination of giant switch statements in drivers
 */

void gfc2d_vi_fpp_set_default(SolverOptions* options) {
  /* VI_FPP doesn't use internal solvers, but we call the VI set_default for proper
   * initialization */
  variationalInequality_FixedPointProjection_set_default(options);
}

static int gfc2d_vi_fpp_init_wrap(void* problem, SolverOptions* options) {
  /* set_default already called by solver_options_create */
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int gfc2d_vi_fpp_solve_wrap(void* problem, double* reaction, double* velocity,
                                   double* globalVelocity, SolverOptions* options) {
  int info = NUMERICS_OK;
  gfc2d_VI_FixedPointProjection((GlobalFrictionContactProblem*)problem, reaction, velocity,
                                globalVelocity, &info, options);
  return info;
}

static void gfc2d_vi_fpp_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER_3VAR(
    GFC2D_VI_FPP, "GFC2D_VI_FPP",
    "Variational Inequality Fixed Point Projection for 2D Global Friction Contact",
    gfc2d_vi_fpp_init_wrap, gfc2d_vi_fpp_solve_wrap, gfc2d_vi_fpp_free_wrap,
    NULL,                     /* error function */
    gfc2d_vi_fpp_set_default, /* set_default */
    1000,                     /* default_max_iter */
    1e-4,                     /* default_tol */
    0 /* is_local_solver */);
