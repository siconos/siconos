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
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for pow, fabs, sqrt, isinf
#include <stdlib.h>  // for calloc, malloc, srand
#include <string.h>  // for NULL, memcpy

#include "FrictionContact_options.h"  // for SICONOS_FRICTION_3D_IP...
#include "NumericsArrays.h"           // for uint_shuffle
#include "NumericsFwd.h"              // for SolverOptions, Rolling...
#include "NumericsMatrix.h"
#include "NumericsVector.h"
#include "SiconosBlas.h"    // for cblas_dnrm2
#include "SolverOptions.h"  // for SolverOptions, SICONOS...
#include "SparseBlockMatrix.h"
#include "naming_conventions.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "rolling_fc_Solvers.h"                       // for RollingComputeErrorPtr
#include "rolling_friction_3d_compute_error.h"        // for rolling_friction_3d_compute_e...
#include "rolling_friction_3d_local_problem_tools.h"  // for rolling_friction_3d_local_pro...
#include "rolling_friction_3d_onecone_nonsmooth_Newton_solvers.h"  // for rolling_friction_3d_projectio...
#include "rolling_friction_3d_projection.h"  // for rolling_friction_3d_projectio...
#include "rolling_friction_3d_short_names.h"
#include "siconos_debug.h"  // for DEBUG_PRINTF, DEBUG_BEGIN
#include "solver_registry.h"

// #define FCLIB_OUTPUT

#ifdef FCLIB_OUTPUT
static int fccounter = -1;
#include "fclib_interface.h"
#endif

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void rolling_friction_3d_nsgs_initialize_local_solver(
    RollingSolverPtr *solve, RollingUpdatePtr *update, RollingFreeSolverNSGSPtr *freeSolver,
    RollingComputeErrorPtr *computeError, RollingFrictionContactProblem *problem,
    RollingFrictionContactProblem *localproblem, SolverOptions *options) {
  SolverOptions *local_opts = options->internalSolvers[0];
  /** Connect to local solver */

  switch (local_opts->solverId) {
    case SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration: {
      *solve = &rolling_friction_3d_projectionOnConeWithLocalIteration_solve;
      *update = &rolling_friction_3d_projection_update;
      *freeSolver =
          (RollingFreeSolverNSGSPtr)&rolling_friction_3d_projectionOnConeWithLocalIteration_free;
      *computeError = &rolling_friction_3d_compute_error;
      rolling_friction_3d_projectionOnConeWithLocalIteration_initialize(problem, localproblem,
                                                                        local_opts);
      break;
    }
    case SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone: {
      *solve = &rolling_friction_3d_projectionOnCone_solve;
      *update = &rolling_friction_3d_projection_update;
      *freeSolver = (RollingFreeSolverNSGSPtr)&rolling_friction_3d_projection_free;
      *computeError = &rolling_friction_3d_compute_error;
      rolling_friction_3d_projection_initialize(problem, localproblem);
      break;
    }
    case SICONOS_ROLLING_FRICTION_3D_ONECONTACT_NSN: {
      *solve = &rolling_friction_3d_onecone_nonsmooth_Newton_solvers_solve;
      *update = &rolling_friction_3d_onecone_nonsmooth_Newton_update;
      *freeSolver =
          (RollingFreeSolverNSGSPtr)&rolling_friction_3d_onecone_nonsmooth_Newton_solvers_free;
      *computeError = &rolling_friction_3d_compute_error;
      rolling_friction_3d_onecone_nonsmooth_Newton_solvers_initialize(
          problem, localproblem, options->internalSolvers[0]);
      break;
    }
    default: {
      numerics_error(
          "rolling_friction_3d_nsgs_initialize_local_solver",
          "Numerics, rolling_friction_3d_nsgs failed. Unknown internal solver : %s.\n",
          solver_options_id_to_name(local_opts->solverId));
    }
  }
}

// Never used ...
/* static */
/* unsigned int* allocShuffledContacts(RollingFrictionContactProblem *problem, */
/*                                     SolverOptions *options) */
/* { */
/*   unsigned int *scontacts = 0; */
/*   unsigned int nc = problem->numberOfContacts; */
/*   if(options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
 * SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE|| */
/*       options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
 * SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP) */
/*   { */
/*     if(options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED] >0) */
/*     { */
/*       srand((unsigned int)options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED]); */
/*     } */
/*     else */
/*       srand(1); */
/*     scontacts = (unsigned int *) malloc(nc * sizeof(unsigned int)); */
/*     for(unsigned int i = 0; i < nc ; ++i) */
/*     { */
/*       scontacts[i] = i; */
/*     } */
/*     uint_shuffle(scontacts, nc); */
/*   } */
/*   return scontacts; */
/* } */

static unsigned int *allocfreezingContacts(RollingFrictionContactProblem *problem,
                                           SolverOptions *options) {
  unsigned int *fcontacts = 0;
  unsigned int nc = problem->numberOfContacts;
  if (options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
    fcontacts = (unsigned int *)malloc(nc * sizeof(unsigned int));
    for (unsigned int i = 0; i < nc; ++i) {
      fcontacts[i] = 0;
    }
  }
  return fcontacts;
}

static int solveLocalReaction(RollingUpdatePtr update_localproblem,
                              RollingSolverPtr local_solver, unsigned int contact,
                              RollingFrictionContactProblem *problem,
                              RollingFrictionContactProblem *localproblem, double *reaction,
                              SolverOptions *local_opts, double r_local[5]) {
  (*update_localproblem)(contact, problem, localproblem, reaction, local_opts);

  local_opts->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER] = contact;

  r_local[0] = reaction[contact * 5 + 0];
  r_local[1] = reaction[contact * 5 + 1];
  r_local[2] = reaction[contact * 5 + 2];
  r_local[3] = reaction[contact * 5 + 3];
  r_local[4] = reaction[contact * 5 + 4];

  return (*local_solver)(localproblem, r_local, local_opts);
}

static void performRelaxation(double r_local[5], double *old_r_local, double omega) {
  r_local[0] = omega * r_local[0] + (1.0 - omega) * old_r_local[0];
  r_local[1] = omega * r_local[1] + (1.0 - omega) * old_r_local[1];
  r_local[2] = omega * r_local[2] + (1.0 - omega) * old_r_local[2];
  r_local[3] = omega * r_local[3] + (1.0 - omega) * old_r_local[3];
  r_local[4] = omega * r_local[4] + (1.0 - omega) * old_r_local[4];
}

static double light_error_squared(double r_local[5], double *old_r_local) {
  return (pow(old_r_local[0] - r_local[0], 2) + pow(old_r_local[1] - r_local[1], 2) +
          pow(old_r_local[2] - r_local[2], 2) + pow(old_r_local[3] - r_local[3], 2) +
          pow(old_r_local[4] - r_local[4], 2));
}

static double squared_norm(double r_local[5]) {
  return (pow(r_local[0], 2) + pow(r_local[1], 2) + pow(r_local[2], 2) + pow(r_local[3], 2) +
          pow(r_local[4], 2));
}

static void acceptLocalReactionFiltered(RollingFrictionContactProblem *localproblem,
                                        SolverOptions *local_opts, unsigned int contact,
                                        unsigned int iter, double *reaction,
                                        double r_local[5]) {
  if (isnan(SOLVER_RESIDUAL(local_opts)) || isinf(SOLVER_RESIDUAL(local_opts)) ||
      SOLVER_RESIDUAL(local_opts) > 1.0) {
    DEBUG_EXPR(rollingFrictionContact_display(localproblem));
    DEBUG_PRINTF(
        "Discard local reaction for contact %i at iteration %i "
        "with local_error = %e\n",
        contact, iter, SOLVER_RESIDUAL(local_opts));

    numerics_printf(
        "Discard local reaction for contact %i at iteration %i "
        "with local_error = %e",
        contact, iter, SOLVER_RESIDUAL(local_opts));
  } else
    memcpy(&reaction[contact * 5], r_local, sizeof(double) * 5);
}

static void acceptLocalReactionUnconditionally(unsigned int contact, double *reaction,
                                               double r_local[5]) {
  memcpy(&reaction[contact * 5], r_local, sizeof(double) * 5);
}

static double calculateLightError(double light_error_sum, unsigned int nc, double *reaction,
                                  double *norm_r) {
  DEBUG_BEGIN("calculateLightError(...)\n");
  double error = sqrt(light_error_sum);
  *norm_r = cblas_dnrm2(nc * 5, reaction, 1);
  if (fabs(*norm_r) > DBL_EPSILON) error /= (*norm_r);
  DEBUG_PRINTF("error = %f\n", error);
  DEBUG_END("calculateLightError(...)\n");
  return error;
}

static double calculateFullErrorAdaptiveInterval(RollingFrictionContactProblem *problem,
                                                 RollingComputeErrorPtr computeError,
                                                 SolverOptions *options, int iter,
                                                 double *reaction, double *velocity,
                                                 double tolerance, double norm_q) {
  double error = 1e+24;
  if (options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] > 0) {
    if (iter % options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] == 0) {
      (*computeError)(problem, reaction, velocity, tolerance, options, norm_q, &error);
      if (error > tolerance && options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                                   SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE)
        options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] *= 2;
    }
    numerics_printf(
        "--------------- RFC3D - NSGS - Iteration %i "
        "options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = %i, "
        "options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] = % i",
        iter, options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY],
        options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION]);
  } else
    (*computeError)(problem, reaction, velocity, tolerance, options, norm_q, &error);

  return error;
}

static double calculateFullErrorFinal(RollingFrictionContactProblem *problem,
                                      SolverOptions *options,
                                      RollingComputeErrorPtr computeError, double *reaction,
                                      double *velocity, double tolerance, double norm_q) {
  double absolute_error;
  (*computeError)(problem, reaction, velocity, tolerance, options, norm_q, &absolute_error);

  if (verbose > 0) {
    if (absolute_error > SOLVER_TOL(options)) {
      numerics_printf(
          "------- RFC3D - NSGS - Warning absolute "
          "Residual = %14.7e is larger than required precision = %14.7e",
          absolute_error, SOLVER_TOL(options));
    } else {
      numerics_printf(
          "------- RFC3D - NSGS - absolute "
          "Residual = %14.7e is smaller than required precision = %14.7e",
          absolute_error, SOLVER_TOL(options));
    }
  }
  return absolute_error;
}

static int determine_convergence(double error, double tolerance, int iter,
                                 SolverOptions *options) {
  int hasNotConverged = 1;
  if (error < tolerance) {
    hasNotConverged = 0;
    numerics_printf("---- RFC3D - NSGS - | %3d | %14.7e | %7.3e |", iter, error, tolerance);
  } else {
    numerics_printf("---- RFC3D - NSGS - | %3d | %14.7e | %7.3e |", iter, error, tolerance);
  }
  return hasNotConverged;
}

static int determine_convergence_with_full_final(RollingFrictionContactProblem *problem,
                                                 SolverOptions *options,
                                                 RollingComputeErrorPtr computeError,
                                                 double *reaction, double *velocity,
                                                 double *tolerance, double norm_q,
                                                 double error, int iter) {
  int hasNotConverged = 1;
  if (error < *tolerance) {
    hasNotConverged = 0;
    numerics_printf("---- RFC3D - NSGS - | %3d | %14.7e | %7.3e |", iter, error, *tolerance);

    double absolute_error = calculateFullErrorFinal(problem, options, computeError, reaction,
                                                    velocity, SOLVER_TOL(options), norm_q);

    if (absolute_error > SOLVER_TOL(options)) {
      if (error < DBL_EPSILON) {
        /* in this case, the relative error is very small
           (meaning that the nsgs loop does not
           improve accuracy).
           We try to tighten the local solver tolerance */
        SET_LOCAL_SOLVER_TOL(
            options->internalSolvers[0],
            fmax(LOCAL_SOLVER_TOL(options->internalSolvers[0]) / 100., DBL_EPSILON * 1e-6));
        numerics_printf(
            "------- RFC3D - NSGS - We modify the local solver tolerance precision to reach "
            "accuracy to %e",
            LOCAL_SOLVER_TOL(options->internalSolvers[0]));

      } else {
        *tolerance = error / absolute_error * SOLVER_TOL(options);
        assert(*tolerance > 0.0 && "tolerance has to be positive");
        numerics_printf(
            "------- RFC3D - NSGS - We modify the required incremental precision to reach "
            "accuracy to %e",
            *tolerance);
      }
      hasNotConverged = 1;
    } else {
      numerics_printf(
          "------- RFC3D - NSGS - The incremental precision is sufficient to reach accuracy "
          "to %e",
          *tolerance);
    }

  } else {
    numerics_printf("---- RFC3D - NSGS - | %3d | %14.7e | %7.3e |", iter, error, *tolerance);
  }
  return hasNotConverged;
}

static void statsIterationCallback(RollingFrictionContactProblem *problem,
                                   SolverOptions *options, double *reaction, double *velocity,
                                   double error) {
  if (options->callback) {
    options->callback->collectStatsIteration(options->callback->env,
                                             problem->numberOfContacts * 5, reaction, velocity,
                                             error, NULL);
  }
}

void rolling_friction_3d_nsgs(RollingFrictionContactProblem *problem, double *reaction,
                              double *velocity, int *info, SolverOptions *options) {
  /* problem->mu_r[0]=0.1; */
  /* problem->mu[0]=1.0; */

  /* for (int c = 0; c < problem->numberOfContacts; c++) { */
  /*   problem->mu[0]=0.0;  */
  /*   }     */

  /* verbose=1; */
  /* int and double parameters */
  int *iparam = options->iparam;
  double *dparam = options->dparam;

  /* Number of contacts */
  unsigned int nc = problem->numberOfContacts;

  /* Maximum number of iterations */
  int itermax = SOLVER_MAX_ITER(options);

  /* Tolerance */
  double tolerance = SOLVER_TOL(options);
  double local_tolerance_save = LOCAL_SOLVER_TOL(options->internalSolvers[0]);
  double norm_q = cblas_dnrm2(nc * 5, problem->q, 1);
  double omega = dparam[SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE];

  double norm_r[] = {1e24};
  RollingSolverPtr local_solver = NULL;
  RollingUpdatePtr update_localproblem = NULL;
  RollingFreeSolverNSGSPtr freeSolver = NULL;
  RollingComputeErrorPtr computeError = NULL;

  RollingFrictionContactProblem *localproblem;
  double r_local[5];

  /*****  NSGS Iterations *****/
  int iter = 0;      /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  unsigned int contact; /* Number of the current row of blocks in M */
  unsigned int *scontacts = NULL;
  unsigned int *freeze_contacts = NULL;

  SparseBlockStructuredMatrix *matrix1 = problem->M->matrix1;
  if (problem->M->storageType == NM_SPARSE) {
    if (problem->M->matrix1) {
      printf("Warning matrix 1 different from NULL");
    }

    problem->M->matrix1 = NM_extract_diagonal_blocks(problem->M, problem->dimension);
  }
  /* Solver initialization continues below */

  if (options->numberOfInternalSolvers < 1) {
    numerics_error("rolling_friction_3d_nsgs",
                   "The NSGS method needs options for the internal solvers, "
                   "options[0].numberOfInternalSolvers should be >= 1");
  }
  SolverOptions *local_opts = options->internalSolvers[0];

  /*****  Initialize various solver options *****/
  localproblem = rolling_friction_3d_local_problem_allocate(problem);

  rolling_friction_3d_nsgs_initialize_local_solver(
      &local_solver, &update_localproblem, (RollingFreeSolverNSGSPtr *)&freeSolver,
      &computeError, problem, localproblem, options);

  freeze_contacts = allocfreezingContacts(problem, options);

  /*****  Check solver options *****/
  if (!(iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE ||
        iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE ||
        iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
            SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP)) {
    numerics_error("rolling_friction_3d_nsgs",
                   "iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] must be equal to "
                   "SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE (0), "
                   "SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE (1) or "
                   "SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP (2)");
    return;
  }

  if (!(iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
            SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL ||
        iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
            SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL ||
        iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
            SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT ||
        iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
            SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE)) {
    numerics_error("rolling_friction_3d_nsgs",
                   "iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] must be equal to "
                   "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL (0), "
                   "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL (1), "
                   "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT (2) or "
                   "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE (3)");
    return;
  }

  /*****  NSGS Iterations *****/

  /* A special case for the most common options (should correspond
   * with mechanics_run.py **/
  if (iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE &&
      iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] == 0 &&
      iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE &&
      iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] ==
          SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE &&
      iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
          SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT) {
    while ((iter < itermax) && (hasNotConverged > 0)) {
      ++iter;
      double light_error_sum = 0.0;

      rolling_friction_3d_set_internalsolver_tolerance(problem, options, local_opts, error);

      for (unsigned int i = 0; i < nc; ++i) {
        contact = i;

        solveLocalReaction(update_localproblem, local_solver, contact, problem, localproblem,
                           reaction, local_opts, r_local);

        light_error_sum += light_error_squared(r_local, &reaction[contact * 5]);

        /* #if 0 */
        acceptLocalReactionFiltered(localproblem, local_opts, contact, iter, reaction,
                                    r_local);
      }

      error = calculateLightError(light_error_sum, nc, reaction, norm_r);

      hasNotConverged = determine_convergence(error, tolerance, iter, options);

      statsIterationCallback(problem, options, reaction, velocity, error);
    }
  }

  /* All other cases, we put all the ifs inline.. otherwise, too many
   * variations to have dedicated loops, but add more if there are
   * common cases to avoid checking booleans on every iteration. **/
  else {
    double *light_error_2 = calloc(nc, sizeof(double));

    while ((iter < itermax) && (hasNotConverged > 0)) {
      ++iter;
      double light_error_sum = 0.0;

      rolling_friction_3d_set_internalsolver_tolerance(problem, options, local_opts, error);
      unsigned int number_of_freezed_contact = 0;

      double tmp_criteria1 = tolerance * tolerance / (nc * nc * 1000);
      double tmp_criteria2 = *norm_r * *norm_r / (nc * nc * 1000);
      if (iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
        for (unsigned int i = 0; i < nc; ++i) {
          if (freeze_contacts[i] > 0) number_of_freezed_contact++;
        }
        if (number_of_freezed_contact >= nc - 1) {
          numerics_printf_verbose(1, "number of freezed contact is too large : %i\n",
                                  number_of_freezed_contact);
          for (unsigned int c = 0; c < nc; ++c) freeze_contacts[c] = 0;
        }
      }
      numerics_printf_verbose(2, "Number of freezed contacts  = %i",
                              number_of_freezed_contact);

      for (unsigned int i = 0; i < nc; ++i) {
        if (iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
                SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE ||
            iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
                SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP) {
          if (iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
              SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP)
            uint_shuffle(scontacts, nc);
          contact = scontacts[i];
        } else
          contact = i;

        if (iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
          if (freeze_contacts[contact] > 0) {
            /* we skip freeze contacts */
            freeze_contacts[contact] -= 1;
            // printf("Contact % i is freezed for %i remaining steps\n", contact,
            // freeze_contacts[contact]);
            light_error_sum += light_error_2[contact];
            continue;
          }
        }

        solveLocalReaction(update_localproblem, local_solver, contact, problem, localproblem,
                           reaction, local_opts, r_local);

        if (iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION] ==
            SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE)
          performRelaxation(r_local, &reaction[contact * 5], omega);

        light_error_2[contact] = light_error_squared(r_local, &reaction[contact * 5]);
        light_error_sum += light_error_2[contact];

        /* int test =100; */
        /* if (contact == test) */
        /* { */
        /*   printf("reaction[%i] = %16.8e\t",3*contact-1,reaction[3*contact]); */
        /*   printf("r_local[%i] = %16.8e\n",2,r_local[0]); */
        /* } */

        if (iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
          double squared_norm_localreaction = squared_norm(r_local);
          int relative_convergence_criteria =
              light_error_2[contact] <= tmp_criteria1 * squared_norm_localreaction;
          int small_reaction_criteria = squared_norm_localreaction <= tmp_criteria2;

          if ((relative_convergence_criteria || small_reaction_criteria) && iter >= 10) {
            /* we  freeze the contact for n iterations*/
            freeze_contacts[contact] =
                options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT];

            // DEBUG_EXPR(
            NV_display(r_local, 5);
            NV_display(&reaction[contact * 5], 5);
            printf("light_error_2 = %e\n", light_error_2[contact]);
            printf("tmp_criteria1 = %e\n", tmp_criteria1);
            printf("tmp_criteria2 = %e\n", tmp_criteria2);
            printf(
                "first criteria relative_convergence_criteria : light_error_2 <= "
                "tmp_criteria1 * squared_norm_localreaction ==> %e <= %e, bool =%i\n",
                light_error_2[contact], tmp_criteria1 * squared_norm_localreaction,
                relative_convergence_criteria);
            printf(
                "second criteria :  squared_norm_localreaction <= tmp_criteria2 ==> %e "
                "<= %e, bool =%i \n",
                squared_norm_localreaction, tmp_criteria2, small_reaction_criteria);
            printf("Contact % i is freezed for %i steps\n", contact,
                   options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT]);
            //);
          }
        }
        if (iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] ==
            SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE)
          acceptLocalReactionFiltered(localproblem, local_opts, contact, iter, reaction,
                                      r_local);
        else
          acceptLocalReactionUnconditionally(contact, reaction, r_local);
      }

      if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
          SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT) {
        error = calculateLightError(light_error_sum, nc, reaction, norm_r);
        hasNotConverged = determine_convergence(error, tolerance, iter, options);
      } else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                 SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
        error = calculateLightError(light_error_sum, nc, reaction, norm_r);
        hasNotConverged =
            determine_convergence_with_full_final(problem, options, computeError, reaction,
                                                  velocity, &tolerance, norm_q, error, iter);

      } else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                 SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL) {
        error = calculateFullErrorAdaptiveInterval(problem, computeError, options, iter,
                                                   reaction, velocity, tolerance, norm_q);
        hasNotConverged = determine_convergence(error, tolerance, iter, options);
      }

      statsIterationCallback(problem, options, reaction, velocity, error);
    }
    /* if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0) */
    /* { */
    /*   int frozen_contact=0; */
    /*   for(unsigned int i = 0 ; i < nc ; ++i) */
    /*   { */
    /*     if (freeze_contacts[contact] >0) */
    /*     { */
    /*       frozen_contact++; */
    /*     } */
    /*   } */
    /*   printf("number of frozen contacts : %i\n", frozen_contact ); */
    /* } */
    free(light_error_2);
  }

  /* Full criterium */
  if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
      SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
    error = calculateFullErrorFinal(problem, options, computeError, reaction, velocity,
                                    tolerance, norm_q);

    hasNotConverged = determine_convergence(error, SOLVER_TOL(options), iter, options);
  }

  *info = hasNotConverged;

  /** return parameter values */
  /* dparam[SICONOS_DPARAM_TOL] = tolerance; */
  SET_SOLVER_RESIDUAL(options, error);
  SET_SOLVER_ITER_DONE(options, iter);
  SET_LOCAL_SOLVER_TOL(options->internalSolvers[0], local_tolerance_save);

  if (problem->M->storageType == NM_SPARSE) {
    SBM_clear_block(problem->M->matrix1);
    SBM_clear(problem->M->matrix1);
    problem->M->matrix1 = matrix1;
  }

  /** Free memory **/
  (*freeSolver)(problem, localproblem, local_opts);
  rolling_friction_3d_local_problem_free(localproblem, problem);
  if (scontacts) free(scontacts);
}

void rolling_friction_3d_nsgs_set_default(SolverOptions *options) {
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] =
      SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] =
      SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE;
  /* options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] =
   * SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE; */
  /* options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] =
   * SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT; */
  options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] = SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE;
  options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED] = 0;
  options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 0;
  options->iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] =
      SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_FALSE;
  options->iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION] =
      SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 0;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
  options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] = 10.0;
  if (options->numberOfInternalSolvers == 0) {
    options->numberOfInternalSolvers = 1;
    options->internalSolvers = calloc(1, sizeof(SolverOptions *));
  }
  assert(options->numberOfInternalSolvers == 1);
  options->internalSolvers[0] =
      solver_options_create(SICONOS_ROLLING_FRICTION_3D_ONECONTACT_NSN);
}

/* Solver registration wrapper functions */
static int rolling_friction_3d_nsgs_init_wrap(void *problem, SolverOptions *options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int rolling_friction_3d_nsgs_solve_wrap(void *problem, double *reaction,
                                               double *velocity, SolverOptions *options) {
  int info = NUMERICS_OK;
  rolling_friction_3d_nsgs((RollingFrictionContactProblem *)problem, reaction, velocity, &info,
                           options);
  return info;
}

static void rolling_friction_3d_nsgs_free_wrap(void *problem, SolverOptions *options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(RFC3D_NSGS, "RFC3D_NSGS",
                "Non-smooth Gauss-Seidel for 3D Rolling Friction Contact",
                rolling_friction_3d_nsgs_init_wrap, rolling_friction_3d_nsgs_solve_wrap,
                rolling_friction_3d_nsgs_free_wrap, NULL,   /* error function */
                rolling_friction_3d_nsgs_set_default, 1000, /* default_max_iter */
                1e-4,                                       /* default_tol */
                0 /* is_local_solver */);
