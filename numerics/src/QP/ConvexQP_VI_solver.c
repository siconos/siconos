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

#include <stdio.h>   // for printf
#include <stdlib.h>  // for free, malloc

#include "ConvexQP.h"                       // for ConvexQP
#include "ConvexQP_Solvers.h"               // for convexQP_VI_solver
#include "ConvexQP_as_VI.h"                 // for ConvexQP_as_VI, Function_...
#include "ConvexQP_computeError.h"          // for convexQP_compute_error_re...
#include "ConvexQP_cst.h"                   // for SICONOS_CONVEXQP_VI_EG
#include "NumericsFwd.h"                    // for SolverOptions, Variationa...
#include "SiconosBlas.h"                    // for cblas_dnrm2
#include "SolverOptions.h"                  // for SolverOptions, SICONOS_DP...
#include "VariationalInequality.h"          // for VariationalInequality
#include "VariationalInequality_Solvers.h"  // for variationalInequality_Ext...
#include "numerics_verbose.h"

/* Solver registration system */
#include "solver_registry.h"
#include "numerics_errors.h"

void convexQP_VI_solver(ConvexQP* problem, double* z, double* w, int* info,
                        SolverOptions* options);

void convexQP_VI_solver(ConvexQP* problem, double* z, double* w, int* info,
                        SolverOptions* options) {
  NumericsMatrix* A = problem->A;
  if (A) {
    *info = numerics_error(
        "ConvexQP_VI_Solver",
        "This solver does not support a specific matrix A different from the identity");
    return;
  }
  /* Dimension of the problem */
  int n = problem->size;

  VariationalInequality* vi = (VariationalInequality*)malloc(sizeof(VariationalInequality));

  // vi.self = &vi;
  vi->F = &Function_VI_CQP;
  vi->ProjectionOnX = &Projection_VI_CQP;

  double error = 1e24;

  ConvexQP_as_VI* convexQP_as_vi = (ConvexQP_as_VI*)malloc(sizeof(ConvexQP_as_VI));
  vi->env = convexQP_as_vi;
  vi->size = n;

  /*set the norm of the VI to the norm of problem->q  */
  double norm_q = cblas_dnrm2(n, problem->q, 1);
  vi->normVI = norm_q;
  vi->istheNormVIset = 1;

  convexQP_as_vi->vi = vi;
  convexQP_as_vi->cqp = problem;

  if (options->solverId == SICONOS_CONVEXQP_VI_FPP) {
    variationalInequality_FixedPointProjection(vi, z, w, info, options);
  } else if (options->solverId == SICONOS_CONVEXQP_VI_EG) {
    variationalInequality_ExtraGradient(vi, z, w, info, options);
  }

  /* **** Criterium convergence **** */
  convexQP_compute_error_reduced(problem, z, w, options->dparam[SICONOS_DPARAM_TOL], options,
                                 norm_q, &error);

  /* for (i =0; i< n ; i++) */
  /* { */
  /*   printf("reaction[%i]=%f\t",i,reaction[i]);
   * printf("velocity[%i]=F[%i]=%f\n",i,i,velocity[i]); */
  /* } */

  if (verbose > 0) {
    if (options->solverId == SICONOS_CONVEXQP_VI_FPP) {
      printf(
          "--------------- CONVEXQP - VI solver (VI_FPP) - #Iteration %i Final Residual = "
          "%14.7e\n",
          options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
    } else if (options->solverId == SICONOS_CONVEXQP_VI_EG) {
      printf(
          "--------------- CONVEXQP - VI solver (VI_EG) - #Iteration %i Final Residual = "
          "%14.7e\n",
          options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
    }
  }
  free(vi);

  free(convexQP_as_vi);
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers SICONOS_CONVEXQP_VI_FPP and SICONOS_CONVEXQP_VI_EG in the
 * global solver registry, enabling:
 * - Dynamic solver lookup by ID
 * - Runtime solver introspection
 * - Elimination of giant switch statements in drivers
 */

static void convexqp_vi_fpp_set_default(SolverOptions* options) {
  variationalInequality_FixedPointProjection_set_default(options);
}

static int convexqp_vi_fpp_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int convexqp_vi_fpp_solve_wrap(void* problem, double* z, double* w, SolverOptions* options) {
  int info = NUMERICS_OK;
  options->solverId = SICONOS_CONVEXQP_VI_FPP;
  convexQP_VI_solver((ConvexQP*)problem, z, w, &info, options);
  return info;
}

static void convexqp_vi_fpp_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_CONVEXQP_VI_FPP, "CONVEXQP_VI_FPP",
                "Variational Inequality Fixed Point Projection for Convex QP",
                convexqp_vi_fpp_init_wrap,
                convexqp_vi_fpp_solve_wrap,
                convexqp_vi_fpp_free_wrap,
                NULL,  /* error function */
                convexqp_vi_fpp_set_default,
                1000,  /* default_max_iter */
                1e-4,  /* default_tol */
                0      /* is_local_solver */);

static void convexqp_vi_eg_set_default(SolverOptions* options) {
  variationalInequality_ExtraGradient_set_default(options);  
}

static int convexqp_vi_eg_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int convexqp_vi_eg_solve_wrap(void* problem, double* z, double* w, SolverOptions* options) {
  int info = NUMERICS_OK;
  options->solverId = SICONOS_CONVEXQP_VI_EG;
  convexQP_VI_solver((ConvexQP*)problem, z, w, &info, options);
  return info;
}

static void convexqp_vi_eg_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_CONVEXQP_VI_EG, "CONVEXQP_VI_EG",
                "Variational Inequality Extra Gradient for Convex QP",
                convexqp_vi_eg_init_wrap,
                convexqp_vi_eg_solve_wrap,
                convexqp_vi_eg_free_wrap,
                NULL,  /* error function */
                convexqp_vi_eg_set_default,
                1000,  /* default_max_iter */
                1e-4,  /* default_tol */
                0      /* is_local_solver */);
