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
#include <assert.h>  // for assert
#include <stdio.h>   // for printf
#include <stdlib.h>  // for malloc, free

#include "FrictionContactProblem.h"        // for FrictionContactProblem
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NumericsFwd.h"                   // for NumericsMatrix
#include "SiconosBlas.h"                   // for cblas_dnrm2
#include "SolverOptions.h"                 // for SICONOS_DPARAM_TOL
#include "fc2d_Solvers.h"                  // for fc2d_nsgs, fc2d_nsgs_set_default
#include "fc3d_short_names.h"
#include "gfc2d_Solvers.h"        // for gfc2d_nsgs_wr
#include "gfc2d_compute_error.h"  // for gfc2d_compute_error
#include "numerics_errors.h"
#include "numerics_verbose.h"
/* Solver registration system */
#include "solver_registry.h"

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include <float.h>
#include <string.h>
#include <time.h>

#include "siconos_debug.h"  // for DEBUG_EXPR, DEBUG_PRINTF

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void gfc2d_nsgs_wr(GlobalFrictionContactProblem* problem, double* reaction, double* velocity,
                   double* globalVelocity, int* info, SolverOptions* options) {
  DEBUG_BEGIN("gfc2d_nsgs_wr\n");
  NumericsMatrix* H = problem->H;
  // We compute only if the local problem has contacts
  DEBUG_PRINTF("Number of contacts = %i \n", H->size1 / problem->dimension);
  if (H->size1 > 0) {
    // Reformulation
    numerics_printf_verbose(1,
                            "Reformulation info a reduced problem onto local variables ...\n");
    FrictionContactProblem* reduced_problem =
        globalFrictionContact_reduce_to_FrictionContact(problem);

    DEBUG_EXPR(frictionContactProblem_display(reduced_problem););
    if (verbose) {
      printf("Call to the fc2d solver ...\n");
    }
    // call nsgs solver for the local problem
    fc2d_nsgs(reduced_problem, reaction, velocity, info, options);

    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);

    /* Number of contacts */
    int nc = problem->numberOfContacts;
    /* Dimension of the problem */
    int m = problem->dimension * nc;
    int n = problem->M->size0;
    double norm_q = cblas_dnrm2(n, problem->q, 1);
    double norm_b = cblas_dnrm2(m, problem->b, 1);
    double error;

    gfc2d_compute_error(problem, reaction, velocity, globalVelocity,
                        options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &error);
    options->dparam[SICONOS_DPARAM_RESIDU] = error;

    frictionContactProblem_free(reduced_problem);
  } else {
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    *info = 0;
  }
  DEBUG_END("gfc2d_nsgs_wr\n");
}

REGISTER_SOLVER_3VAR(GFC2D_NSGS_WR, "GFC2D_NSGS_WR",
                     "Non-smooth Gauss-Seidel for Global 2D Friction Contact with reduction",
                     NULL, NULL, NULL, NULL, /* error function */
                     fc2d_nsgs_set_default,  /* set_default */
                     1000,                   /* default_max_iter */
                     1e-4,                   /* default_tol */
                     0 /* is_local_solver */)
