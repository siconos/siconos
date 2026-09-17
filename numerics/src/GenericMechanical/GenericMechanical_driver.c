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

/*! \file GenericMechanical_driver.c
 *  \brief Driver and default-option setup for GenericMechanicalProblem.
 */

#include <assert.h>  // for assert
#include <stdio.h>   // for size_t
#include <stdlib.h>  // for calloc
#include <string.h>  // for strcat, strcpy

#include "FrictionContact_options.h"    // for SICONOS_FRICTION_3D_ONECON...
#include "GMPReduced.h"                 // for gmp_as_mlcp, gmp_reduced_e...
#include "GenericMechanicalProblem.h"   // for GenericMechanicalProblem
#include "GenericMechanical_Solvers.h"  // for gmp_compute_error, gmp_driver
#include "GenericMechanical_cst.h"      // for SICONOS_GENERIC_MECHANICAL...
#include "Relay_options.h"              // for SICONOS_RELAY_LEMKE
#include "SolverOptions.h"              // for SolverOptions, solver_opti...
#include "lcp_cst.h"                    // for SICONOS_LCP_LEMKE
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "solver_registry.h"

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"  // for DEBUG_BEGIN, DEBUG_END, DEBUG_EXPR

/* See allowed values for iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED]
   in GenericalMechanical_cst.h (GENERIC_MECHANICAL_ISREDUCED enum)
*/
int gmp_driver(GenericMechanicalProblem* problem, double* reaction, double* velocity,
               SolverOptions* options) {
  DEBUG_BEGIN("gmp_driver(...)\n");
  int info = 0;
  DEBUG_EXPR(
      // NM_display(problem->M);
      genericMechanicalProblem_display(problem););
  if (verbose) solver_options_print(options);

  switch (options->solverId) {
    /* Non Smooth Gauss Seidel (NSGS) */
    case SICONOS_GENERIC_MECHANICAL_NSGS: {
      if (options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] ==
          SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS) {
        numerics_printf("gmp_driver : call of gmp_gauss_seidel\n");
        gmp_gauss_seidel(problem, reaction, velocity, &info, options);
      } else if (options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] ==
                 SICONOS_GENERIC_MECHANICAL_SUBS_EQUALITIES) {
        numerics_printf("gmp_driver : call of gmp_reduced_solve\n");
        gmp_reduced_solve(problem, reaction, velocity, &info, options);
      } else if (options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] ==
                 SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES) {
        numerics_printf("gmp_driver : call of gmp_reduced_equality_solve\n");
        gmp_reduced_equality_solve(problem, reaction, velocity, &info, options);
      } else if (options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] ==
                 SICONOS_GENERIC_MECHANICAL_MLCP_LIKE) {
        numerics_printf("gmp_driver : call of mlcp\n");
        gmp_as_mlcp(problem, reaction, velocity, &info, options);
      } else {
        numerics_printf(
            "gmp_driver error, options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED] "
            "wrong value.\n");
      }
      break;
    }
    default: {
      char msg[200];
      strcpy(msg, "Unknown solver : ");
      strcat(msg, solver_options_id_to_name(options->solverId));
      strcat(msg, "\n");
      numerics_warning("gmp_driver", msg);
      return numerics_error("gmp_driver", msg);
    }
  }
  DEBUG_END("gmp_driver(...)\n");
  return info;
}

void gmp_set_default(SolverOptions* options) {
  DEBUG_BEGIN("gmp_set_default(SolverOptions* options)\n");
  /*with Line search 1 without 0.*/
  options->iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_WITH_LINESEARCH] = 0;

  options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
  /*Useful parameter for LS*/
  options->dparam[SICONOS_DPARAM_GMP_COEFF_LS] = 1.0;

  if (options->numberOfInternalSolvers == 0) {
    options->numberOfInternalSolvers = 4;
    options->internalSolvers = calloc(4, sizeof(SolverOptions*));
  } else {
    for (size_t i = 0; i < 4; ++i) solver_options_delete(options->internalSolvers[i]);
  }
  assert(options->numberOfInternalSolvers == 4);

  options->internalSolvers[0] = solver_options_create(SICONOS_LCP_LEMKE);
  options->internalSolvers[1] =
      solver_options_create(SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID);
  options->internalSolvers[2] = solver_options_create(SICONOS_RELAY_LEMKE);
  options->internalSolvers[3] = solver_options_create(SICONOS_FRICTION_2D_NSGS);

  DEBUG_END("gmp_set_default(SolverOptions* options)\n");
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers SICONOS_GENERIC_MECHANICAL_NSGS in the global solver registry, enabling:
 * - Dynamic solver lookup by ID
 * - Runtime solver introspection
 * - Elimination of giant switch statements in drivers
 */

static int gmp_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  /* set_default already called by solver_options_create */
  return NUMERICS_OK;
}

static int gmp_solve_wrap(void* problem, double* reaction, double* velocity,
                          SolverOptions* options) {
  int info = NUMERICS_OK;
  gmp_driver((GenericMechanicalProblem*)problem, reaction, velocity, options);
  return info;
}

static void gmp_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_GENERIC_MECHANICAL_NSGS, "GMP_NSGS",
                "Non-smooth Gauss-Seidel for Generic Mechanical Problem", gmp_init_wrap,
                gmp_solve_wrap, gmp_free_wrap, NULL, /* error function */
                gmp_set_default,                     /* set_default */
                1000,                                /* default_max_iter */
                1e-4,                                /* default_tol */
                0 /* is_local_solver */);
