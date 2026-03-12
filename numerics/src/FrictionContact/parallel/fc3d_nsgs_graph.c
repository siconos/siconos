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
#include <math.h>    // for fabs, sqrt
#include <stdio.h>   // for fclose, fopen
#include <stdlib.h>  // for calloc, malloc
#include <string.h>  // for NULL, memcpy

#include "FrictionContactProblem.h"                    // for FrictionContac...
#include "FrictionContact_options.h"                   // for SICONOS_FRICTI...
#include "Friction_tools.h"                            // for ComputeErrorPtr
#include "fc3d_Solvers.h"                              // for fc3d_nsgs_set_default
#include "NumericsArrays.h"                            // for uint_shuffle
#include "NumericsFwd.h"                               // for SolverOptions
#include "SiconosBlas.h"                               // for cblas_dnrm2
#include "SolverOptions.h"                             // for SolverOptions
#include "fc3d_2NCP_Glocker.h"                         // for NCPGlocker_update
#include "fc3d_NCPGlockerFixedPoint.h"                 // for fc3d_FixedP_in...
#include "fc3d_Path.h"                                 // for fc3d_Path_init...
#include "fc3d_compute_error.h"                        // for fc3d_compute_e...
#include "fc3d_local_problem_tools.h"                  // for fc3d_local_pro...
#include "fc3d_onecontact_nonsmooth_Newton_solvers.h"  // for fc3d_onecontac...
#include "fc3d_projection.h"                           // for fc3d_projectio...
#include "fc3d_unitary_enumerative.h"                  // for fc3d_unitary_e...
#include "fc3d_short_names.h"                  // Short names for solver IDs
#include "numerics_verbose.h"                          // for numerics_printf

/* Solver registration system */
#include "solver_registry.h"
#include "numerics_errors.h"

/* New utility headers for standardized error computation, tolerance management, and naming conventions */
#include "error_computation.h"
#include "tolerance_manager.h"

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */

#include <time.h>

#include "NumericsVector.h"
#include "gfc3d_ipm.h"
#include "op3x3.h"
#include "projectionOnCone.h"  // for projectionOnCone
#include "siconos_debug.h"     // for DEBUG_EXPR

#include "graph_tools.h"

// #define FCLIB_OUTPUT

#ifdef FCLIB_OUTPUT
static int fccounter = -1;
#include "fclib_interface.h"
#endif

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

/* static void fake_compute_error_nsgs(FrictionContactProblem* problem, double *reaction,
 * double *velocity, double tolerance, SolverOptions  *options,  double* error) */
/* { */
/*   int n = 3 * problem->numberOfContacts; */
/*   *error = 0.; */
/*   int i, m; */
/*   m = 5 * n / 3; */
/*   double err = INFINITY; */
/*   for (i = 0 ; i < m ; ++i) */
/*   { */
/*     *error += Compute_NCP_error1(i, err); */
/*   } */
/* } */

static inline void performRelaxation_3(double localreaction[3], double* oldreaction,
                                       double omega) {
  localreaction[0] = omega * localreaction[0] + (1.0 - omega) * oldreaction[0];
  localreaction[1] = omega * localreaction[1] + (1.0 - omega) * oldreaction[1];
  localreaction[2] = omega * localreaction[2] + (1.0 - omega) * oldreaction[2];
}

static inline double light_error_squared_3(double localreaction[3], double* oldreaction) {
  return (pow(oldreaction[0] - localreaction[0], 2) +
          pow(oldreaction[1] - localreaction[1], 2) +
          pow(oldreaction[2] - localreaction[2], 2));
}

static inline double squared_norm_3(double localreaction[3]) {
  return (pow(localreaction[0], 2) + pow(localreaction[1], 2) + pow(localreaction[2], 2));
}

static void fc3d_nsgs_update_parallel(int contact, FrictionContactProblem* problem,
                                      FrictionContactProblem* localproblem, double* reaction,
                                      SolverOptions* options) {
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */
  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  fc3d_local_problem_fill_M(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products
     MLocal.reactionBlock, excluding the block corresponding to the current contact. ****/
  fc3d_local_problem_compute_q_parallel(problem, localproblem, reaction, contact);

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];
}

static void fc3d_nsgs_initialize_local_solver_parallel(
    struct LocalProblemFunctionToolkit* local_function_toolkit, ComputeErrorPtr* computeError,
    FrictionContactProblem* problem, FrictionContactProblem* localproblem,
    SolverOptions* options, SolverOptions* local_opts) {

  *computeError = (ComputeErrorPtr)&fc3d_compute_error;

  if (problem->dimension == 3) {
    local_function_toolkit->copy_local_reaction = cpy3;
    local_function_toolkit->perform_relaxation = &performRelaxation_3;
    local_function_toolkit->light_error_squared = &light_error_squared_3;
    local_function_toolkit->squared_norm = &squared_norm_3;
  }

  /** Connect to local solver */
  switch (local_opts->solverId) {
    /* Projection */
    case OC_PROJ_DIAG: {
      local_function_toolkit->local_solver = &fc3d_projectionWithDiagonalization_solve;
      local_function_toolkit->update_local_problem =
          &fc3d_projectionWithDiagonalization_update;
      local_function_toolkit->free_local_solver = &fc3d_projection_free;
      fc3d_projection_initialize(problem);
      break;
    }
    case OC_PROJ: {
      local_function_toolkit->local_solver = &fc3d_projectionOnCone_solve;
      local_function_toolkit->update_local_problem = &fc3d_nsgs_update_parallel;
      local_function_toolkit->free_local_solver = &fc3d_projection_free;
      fc3d_projection_initialize(problem);
      break;
    }
    case OC_PROJ_LI: {
      local_function_toolkit->local_solver = &fc3d_projectionOnConeWithLocalIteration_solve;
      local_function_toolkit->update_local_problem = &fc3d_nsgs_update_parallel;
      local_function_toolkit->free_local_solver =
          &fc3d_projectionOnConeWithLocalIteration_free;
      fc3d_projectionOnConeWithLocalIteration_initialize(problem,
                                                         local_opts);
      break;
    }
    case OC_PROJ_REG: {
      local_function_toolkit->local_solver = &fc3d_projectionOnCone_solve;
      local_function_toolkit->update_local_problem =
          &fc3d_projection_update_with_regularization_parallel;
      local_function_toolkit->free_local_solver = &fc3d_projection_with_regularization_free;
      fc3d_projection_initialize_with_regularization(problem, localproblem);
      break;
    }
    /* Newton solver (Alart-Curnier) */
    case OC_NSN: {
      local_function_toolkit->local_solver = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      local_function_toolkit->update_local_problem =
          &fc3d_onecontact_nonsmooth_Newton_AC_update_parallel;
      local_function_toolkit->free_local_solver =
          &fc3d_onecontact_nonsmooth_Newton_solvers_free;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem,
                                                          local_opts);
      break;
    }
    case OC_NSN_GP: {
      local_function_toolkit->local_solver = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      local_function_toolkit->update_local_problem =
          &fc3d_onecontact_nonsmooth_Newton_AC_update_parallel;
      local_function_toolkit->free_local_solver =
          &fc3d_onecontact_nonsmooth_Newton_solvers_free;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem,
                                                          local_opts);
      break;
    }
    case OC_NSN_GP_HYBRID: {
      local_function_toolkit->local_solver = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      local_function_toolkit->update_local_problem =
          &fc3d_onecontact_nonsmooth_Newton_AC_update_parallel;
      local_function_toolkit->free_local_solver =
          &fc3d_onecontact_nonsmooth_Newton_solvers_free;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, 
                                                          local_opts);
      break;
    } /* Newton solver (Glocker-Fischer-Burmeister)*/
    case SICONOS_FRICTION_3D_NCPGlockerFBNewton: {
      local_function_toolkit->local_solver = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      local_function_toolkit->update_local_problem = &NCPGlocker_update;
      local_function_toolkit->free_local_solver =
          &fc3d_onecontact_nonsmooth_Newton_solvers_free;
      // *computeError = &fake_compute_error;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, 
                                                          local_opts);
      break;
    }
    /* Path solver (Glocker Formulation) */
    case SICONOS_FRICTION_3D_NCPGlockerFBPATH: {
      local_function_toolkit->local_solver = &fc3d_Path_solve;
      local_function_toolkit->free_local_solver = (FreeLocalSolverPtr)&fc3d_Path_free;
      local_function_toolkit->update_local_problem = &NCPGlocker_update;

      // *computeError = &fake_compute_error;
      fc3d_Path_initialize(problem, localproblem, local_opts);
      break;
    }

    /* Fixed Point solver (Glocker Formulation) */
    case FC3D_NCPG_FP: {
      local_function_toolkit->local_solver = &fc3d_FixedP_solve;
      local_function_toolkit->update_local_problem = &NCPGlocker_update;
      local_function_toolkit->free_local_solver = (FreeLocalSolverPtr)&fc3d_FixedP_free;
      /* *computeError = &fake_compute_error_nsgs; */
      fc3d_FixedP_initialize(problem, localproblem, local_opts);
      break;
    }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder: {
      local_function_toolkit->local_solver = &fc3d_projectionOnCylinder_solve;
      local_function_toolkit->update_local_problem = &fc3d_projectionOnCylinder_update_parallel;
      local_function_toolkit->free_local_solver =
          (FreeLocalSolverPtr)&fc3d_projectionOnCylinder_free;
      *computeError = (ComputeErrorPtr)&fc3d_Tresca_compute_error;
      fc3d_projectionOnCylinder_initialize(problem, localproblem, options);
      break;
    }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration: {
      local_function_toolkit->local_solver =
          &fc3d_projectionOnCylinderWithLocalIteration_solve;
      local_function_toolkit->update_local_problem = &fc3d_projectionOnCylinder_update_parallel;
      local_function_toolkit->free_local_solver =
          (FreeLocalSolverPtr)&fc3d_projectionOnCylinderWithLocalIteration_free;
      *computeError = (ComputeErrorPtr)&fc3d_Tresca_compute_error;
      fc3d_projectionOnCylinderWithLocalIteration_initialize(problem, localproblem, options,
                                                             local_opts);
      break;
    }
    case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC: {
      local_function_toolkit->local_solver = &fc3d_unitary_enumerative_solve;
      local_function_toolkit->update_local_problem = &fc3d_nsgs_update_parallel;
      local_function_toolkit->free_local_solver =
          (FreeLocalSolverPtr)&fc3d_unitary_enumerative_free;
      fc3d_unitary_enumerative_initialize(localproblem);
      break;
    }
    case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU: {
      local_function_toolkit->local_solver = &fc3d_unitary_enumerative_solve;
      local_function_toolkit->update_local_problem = &fc3d_nsgs_update_parallel;
      local_function_toolkit->free_local_solver =
          (FreeLocalSolverPtr)&fc3d_unitary_enumerative_free;
      fc3d_unitary_enumerative_initialize(localproblem);
      break;
    }
    default: {
      numerics_error("fc3d_nsgs_initialize_local_solver",
                     "Numerics, fc3d_nsgs failed. Unknown internal solver : %s.\n",
                     solver_options_id_to_name(local_opts->solverId));
    }
  }
}

static unsigned int* allocShuffledContacts(FrictionContactProblem* problem,
                                           SolverOptions* options) {
  unsigned int* scontacts = 0;
  unsigned int nc = problem->numberOfContacts;
  if (options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
          SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE ||
      options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
          SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP) {
    if (options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED] > 0) {
      srand((unsigned int)options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED]);
    } else
      srand(1);
    scontacts = (unsigned int*)malloc(nc * sizeof(unsigned int));
    for (unsigned int i = 0; i < nc; ++i) {
      scontacts[i] = i;
    }
    uint_shuffle(scontacts, nc);
  }
  return scontacts;
}
static unsigned int* allocfreezingContacts(FrictionContactProblem* problem,
                                           SolverOptions* options) {
  unsigned int* fcontacts = 0;
  unsigned int nc = problem->numberOfContacts;
  if (options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
    fcontacts = (unsigned int*)malloc(nc * sizeof(unsigned int));
    for (unsigned int i = 0; i < nc; ++i) {
      fcontacts[i] = 0;
    }
  }
  return fcontacts;
}

static int solveLocalReaction(UpdatePtr update_localproblem, SolverPtr local_solver,
                              CopyLocalReactionPtr copyLocalReaction, unsigned int contact,
                              FrictionContactProblem* problem,
                              FrictionContactProblem* localproblem, double* reaction,
                              SolverOptions* local_opts, double localreaction[3]) {
  (*update_localproblem)(contact, problem, localproblem, reaction, local_opts);

  local_opts->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER] = contact;

  copyLocalReaction(&(reaction[contact * problem->dimension]), localreaction);

  return (*local_solver)(localproblem, localreaction, local_opts);
}

static void performRelaxation(double localreaction[3], double* oldreaction, double omega) {
  localreaction[0] = omega * localreaction[0] + (1.0 - omega) * oldreaction[0];
  localreaction[1] = omega * localreaction[1] + (1.0 - omega) * oldreaction[1];
  localreaction[2] = omega * localreaction[2] + (1.0 - omega) * oldreaction[2];
}

static double light_error_squared(double localreaction[3], double* oldreaction) {
  return (pow(oldreaction[0] - localreaction[0], 2) +
          pow(oldreaction[1] - localreaction[1], 2) +
          pow(oldreaction[2] - localreaction[2], 2));
}

static double squared_norm(double localreaction[3]) {
  return (pow(localreaction[0], 2) + pow(localreaction[1], 2) + pow(localreaction[2], 2));
}

static int file_exists(const char* fname) {
  FILE* file;
  if ((file = fopen(fname, "r"))) {
    fclose(file);
    return 1;
  }
  return 0;
}

static void acceptLocalReactionFiltered(FrictionContactProblem* localproblem,
                                        SolverOptions* localsolver_options,
                                        unsigned int contact, unsigned int iter,
                                        double* reaction, double localreaction[3]) {
  if (isnan(localsolver_options->dparam[SICONOS_DPARAM_RESIDU]) ||
      isinf(localsolver_options->dparam[SICONOS_DPARAM_RESIDU]) ||
      localsolver_options->dparam[SICONOS_DPARAM_RESIDU] > 1.0) {
    DEBUG_EXPR(frictionContact_display(localproblem));
    DEBUG_PRINTF(
        "Discard local reaction for contact %i at iteration %i "
        "with local_error = %e\n",
        contact, iter, localsolver_options->dparam[SICONOS_DPARAM_RESIDU]);

#ifdef FCLIB_OUTPUT

    /* printf("step counter value = %i\n", localsolver_options->iparam[19]); */
    char fname[256];
    fccounter++;
    snprintf(fname, sizeof(fname), "./local_problem/localproblem_%i_%i.hdf5", contact,
             localsolver_options->iparam[19]);

    if (file_exists(fname)) {
      /* printf(" %s already dumped\n", fname); */
    } else {
      printf("Dump %s\n", fname);
      int n = 100;
      char* title = (char*)malloc(n * sizeof(char));
      strcpy(title, "Bad local problem dump in hdf5");
      char* description = (char*)malloc(n * sizeof(char));
      strcpy(description, "Rewriting in hdf5 from siconos ");
      strcat(description, fname);
      strcat(description, " in FCLIB format");
      char* mathInfo = (char*)malloc(n * sizeof(char));
      strcpy(mathInfo, "unknown");

      frictionContact_fclib_write(localproblem, title, description, mathInfo, fname, 3);

      printf("end of dump %s\n", fname);
      free(title);
      free(description);
      free(mathInfo);
    }

#endif

    numerics_printf(
        "Discard local reaction for contact %i at iteration %i "
        "with local_error = %e",
        contact, iter, localsolver_options->dparam[SICONOS_DPARAM_RESIDU]);
  } else
    memcpy(&reaction[contact * 3], localreaction, sizeof(double) * 3);
}

static void acceptLocalReactionUnconditionally(unsigned int contact, double* reaction,
                                               double localreaction[3]) {
  memcpy(&reaction[contact * 3], localreaction, sizeof(double) * 3);
}

static double calculateLightError(double light_error_sum, unsigned int nc, double* reaction,
                                  double* norm_r) {
  double error = sqrt(light_error_sum);
  *norm_r = cblas_dnrm2(nc * 3, reaction, 1);
  // printf("error = %e, norm_r = %e\n", error, *norm_r);
  if (fabs(*norm_r) > DBL_EPSILON) error /= (*norm_r);
  return error;
}

static double calculateFullErrorAdaptiveInterval(FrictionContactProblem* problem,
                                                 ComputeErrorPtr computeError,
                                                 SolverOptions* options, int iter,
                                                 double* reaction, double* velocity,
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
        "--------------- FC3D - NSGS - Iteration %i "
        "options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = %i, "
        "options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] = % i",
        iter, options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY],
        options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION]);
  } else
    (*computeError)(problem, reaction, velocity, tolerance, options, norm_q, &error);

  return error;
}

static double calculateFullErrorFinal(FrictionContactProblem* problem, SolverOptions* options,
                                      ComputeErrorPtr computeError, double* reaction,
                                      double* velocity, double tolerance, double norm_q) {
  double absolute_error;
  (*computeError)(problem, reaction, velocity, tolerance, options, norm_q, &absolute_error);

  if (verbose > 0) {
    if (absolute_error > options->dparam[SICONOS_DPARAM_TOL]) {
      numerics_printf(
          "------- FC3D - NSGS - Warning absolute "
          "Residual = %14.7e is larger than required precision = %14.7e",
          absolute_error, options->dparam[SICONOS_DPARAM_TOL]);
    } else {
      numerics_printf(
          "------- FC3D - NSGS - absolute "
          "Residual = %14.7e is smaller than required precision = %14.7e",
          absolute_error, options->dparam[SICONOS_DPARAM_TOL]);
    }
  }
  return absolute_error;
}

static int determine_convergence(double error, double tolerance, int iter,
                                 SolverOptions* options) {
  int hasNotConverged = 1;
  if (error < tolerance) {
    hasNotConverged = 0;
    numerics_printf(
        "--------------- FC3D - NSGS - Iteration %i "
        "Residual = %14.7e < %7.3e\n",
        iter, error, tolerance);
  } else {
    numerics_printf(
        "--------------- FC3D - NSGS - Iteration %i "
        "Residual = %14.7e > %7.3e\n",
        iter, error, tolerance);
  }
  return hasNotConverged;
}

static int determine_convergence_with_full_final(FrictionContactProblem* problem,
                                                 SolverOptions* options,
                                                 ComputeErrorPtr computeError,
                                                 double* reaction, double* velocity,
                                                 double* tolerance, double norm_q,
                                                 double error, int iter) {
  int hasNotConverged = 1;
  if (error < *tolerance) {
    hasNotConverged = 0;
    numerics_printf(
        "--------------- FC3D - NSGS - Iteration %i "
        "Residual = %14.7e < %7.3e",
        iter, error, *tolerance);

    double absolute_error =
        calculateFullErrorFinal(problem, options, computeError, reaction, velocity,
                                options->dparam[SICONOS_DPARAM_TOL], norm_q);
    if (absolute_error > options->dparam[SICONOS_DPARAM_TOL]) {
      *tolerance = error / absolute_error * options->dparam[SICONOS_DPARAM_TOL];
      /* assert(*tolerance > 0.0 && "tolerance has to be positive"); */

      /* if (*tolerance < DBL_EPSILON) */
      /* { */
      /*   numerics_warning("determine_convergence_with_full_fina", "We try to set a very smal
       * tolerance"); */
      /*   *tolerance = DBL_EPSILON; */
      /* } */
      numerics_printf(
          "------- FC3D - NSGS - We modify the required incremental precision to reach "
          "accuracy to %e",
          *tolerance);
      hasNotConverged = 1;
    } else {
      numerics_printf(
          "------- FC3D - NSGS - The incremental precision is sufficient to reach accuracy to "
          "%e",
          *tolerance);
    }

  } else {
    numerics_printf(
        "--------------- FC3D - NSGS - Iteration %i "
        "Residual = %14.7e > %7.3e",
        iter, error, *tolerance);
  }
  return hasNotConverged;
}

static void statsIterationCallback(FrictionContactProblem* problem, SolverOptions* options,
                                   double* reaction, double* velocity, double error) {
  if (options->callback) {
    options->callback->collectStatsIteration(options->callback->env,
                                             problem->numberOfContacts * 3, reaction, velocity,
                                             error, NULL);
  }
}

void fc3d_nsgs_graph(FrictionContactProblem* problem, double* reaction, double* velocity,
                     int* info, SolverOptions* options) {
  /* verbose=1; */

  /* Number of contacts */
  unsigned int nc = problem->numberOfContacts;

  /* Maximum number of iterations */
  int itermax = SOLVER_MAX_ITER(options);

  /* Tolerance */
  double norm_q = cblas_dnrm2(nc * 3, problem->q, 1);
  double omega = options->dparam[SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE];

  double norm_r[] = {1e24};
  if (options->numberOfInternalSolvers < 1) {
    numerics_error("fc3d_nsgs",
                   "The NSGS method needs options for the internal solvers, "
                   "options[0].numberOfInternalSolvers should be >= 1");
  }

  /* Get local solver options - use consistent naming */
  SolverOptions* local_opts = options->internalSolvers[0];

  /* Initialize tolerance manager for unified tolerance handling */
  ToleranceManager tol_manager;
  tolerance_manager_init(&tol_manager, SOLVER_TOL(options), local_opts);
  
  /* Working tolerance (may be adapted during iterations) */
  double tolerance = tol_manager.working_tolerance;

  ComputeErrorPtr computeError = NULL;
  
  struct LocalProblemFunctionToolkit* localProblemFunctionToolkit =
      localProblemFunctionToolkit_new();

  FrictionContactProblem* localproblem;

  double localreaction[3];

  /*****  NSGS Iterations *****/
  int iter = 0;      /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  unsigned int contact; /* Number of the current row of blocks in M */
  unsigned int* scontacts = NULL;
  unsigned int* freeze_contacts = NULL;

  SparseBlockStructuredMatrix* matrix1 = problem->M->matrix1;
  if (problem->M->storageType == NM_SPARSE) {
    if (problem->M->matrix1) {
      printf("Warning matrix 1 different from NULL");
    }

    problem->M->matrix1 = NM_extract_diagonal_blocks(problem->M, problem->dimension);
  }

  /*****  Initialize various solver options *****/
  localproblem = fc3d_local_problem_allocate(problem);

  fc3d_nsgs_initialize_local_solver_parallel(localProblemFunctionToolkit, &computeError, problem,
                                             localproblem, options, local_opts);

  scontacts = allocShuffledContacts(problem, options);
  freeze_contacts = allocfreezingContacts(problem, options);

  /*****  Check solver options *****/
  if (!(options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE ||
        options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE ||
        options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
            SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP)) {
    numerics_error("fc3d_nsgs",
                   "iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] must be equal to "
                   "SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE (0), "
                   "SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE (1) or "
                   "SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP (2)");
    return;
  }

  if (!(options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
            SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL ||
        options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
            SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL ||
        options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
            SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT ||
        options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
            SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE)) {
    numerics_error("fc3d_nsgs",
                   "iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] must be equal to "
                   "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL (0), "
                   "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL (1), "
                   "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT (2) or "
                   "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE (3)");
    return;
  }

  /* Coloring */
  size_t n_colors = 0;
  size_t* partition_size = NULL;
  size_t** partitions = NULL;
  color_graph_block(problem->numberOfContacts, problem->M, &n_colors, &partition_size,
                    &partitions);

  /*****  NSGS Iterations *****/

  /* A special case for the most common options (should correspond
   * with mechanics_run.py **/
  if (options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE &&
      options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] == 0 &&
      options->iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION] ==
          SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE &&
      options->iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] ==
          SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE &&
      options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
          SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT) {
    double light_error_sum = 0.;
#pragma omp parallel default(none)                                                        \
    private(contact, localproblem, localreaction, local_opts)                    \
    shared(problem, localProblemFunctionToolkit, computeError, options, \
               iter, itermax, hasNotConverged, error, n_colors, light_error_sum,          \
               partition_size, partitions, reaction, nc, norm_r, tolerance, velocity)
    {
      /* Allocate localproblem for each thread */
      localproblem = fc3d_local_problem_allocate(problem);

      local_opts = solver_options_create(options->internalSolvers[0]->solverId);

      fc3d_nsgs_initialize_local_solver_parallel(localProblemFunctionToolkit, &computeError, problem,
                                                 localproblem, options, local_opts);

      while ((iter < itermax) && (hasNotConverged > 0)) {
        fc3d_set_internalsolver_tolerance(problem, options, local_opts, error);

        for (size_t color = 0; color < n_colors; color++) {
#pragma omp for schedule(static) reduction(+ : light_error_sum)
          for (int v = 0; v < partition_size[color]; v++) {
            contact = partitions[color][v];

            solveLocalReaction(localProblemFunctionToolkit->update_local_problem,
                               localProblemFunctionToolkit->local_solver,
                               localProblemFunctionToolkit->copy_local_reaction, contact, problem,
                               localproblem, reaction, local_opts, localreaction);

            light_error_sum += localProblemFunctionToolkit->light_error_squared(localreaction, &reaction[contact * 3]);

            /* #if 0 */
            acceptLocalReactionFiltered(localproblem, local_opts, contact, iter,
                                        reaction, localreaction);
          }
        }

#pragma omp single
        {
          error = calculateLightError(light_error_sum, nc, reaction, norm_r);

          hasNotConverged = determine_convergence(error, tolerance, iter, options);

          statsIterationCallback(problem, options, reaction, velocity, error);
          ++iter;
          light_error_sum = 0.;
        }
      }
      /* Free stuff allocated at the beginning of parallel ?
      Or don't care because it's done automatically since the parallel region finishes right
      after ?
      */
    }
  }

  /* All other cases, we put all the ifs inline.. otherwise, too many
   * variations to have dedicated loops, but add more if there are
   * common cases to avoid checking booleans on every iteration. **/
  else {
    /* verbose=1; */
    double light_error_sum = 0.;
    double light_error_2;
    double tmp_criteria1, tmp_criteria2;
    int number_of_freezed_contact;

#pragma omp parallel default(none)                                                          \
    private(contact, localproblem, localreaction, light_error_2, local_opts)       \
    shared(problem, localProblemFunctionToolkit, computeError, options,   \
               iter, itermax, hasNotConverged, error, number_of_freezed_contact,            \
               tmp_criteria1, tmp_criteria2, norm_r, nc, n_colors, light_error_sum, \
               partition_size, partitions, reaction, freeze_contacts, scontacts, norm_q,    \
               omega, tolerance, velocity)
    {
      /* Allocate localproblem for each thread */
      localproblem = fc3d_local_problem_allocate(problem);

      local_opts = solver_options_create(options->internalSolvers[0]->solverId);

      fc3d_nsgs_initialize_local_solver_parallel(localProblemFunctionToolkit, &computeError, problem,
                                                 localproblem, options, local_opts);

      while ((iter < itermax) && (hasNotConverged > 0)) {
        fc3d_set_internalsolver_tolerance(problem, options, local_opts, error);

#pragma omp single
        {
          number_of_freezed_contact = 0;
          tmp_criteria1 = tolerance * tolerance * 100 * 100;
          tmp_criteria2 = *norm_r * *norm_r / (nc * nc * 1000);

          if (options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
            for (unsigned int i = 0; i < nc; ++i) {
              if (freeze_contacts[i] > 0) number_of_freezed_contact++;
            }
            if (number_of_freezed_contact >= nc - 1) {
              // printf("number of freezed contact too large\n");
              for (unsigned int c = 0; c < nc; ++c) freeze_contacts[c] = 0;
            }
          }
        }

        for (size_t color = 0; color < n_colors; color++) {
#pragma omp for schedule(static) reduction(+ : light_error_sum)
          for (int v = 0; v < partition_size[color]; v++) {
            int i = partitions[color][v];

            if (options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
                    SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE ||
                options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
                    SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP) {
              if (options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
                  SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP) {
                uint_shuffle(scontacts, nc);
              }
              contact = scontacts[i];
            } else
              contact = i;

            if (options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
              if (freeze_contacts[contact] > 0) {
                /* we skip freeze contacts */
                freeze_contacts[contact] -= 1;
                continue;
              }
            }

            solveLocalReaction(localProblemFunctionToolkit->update_local_problem,
                               localProblemFunctionToolkit->local_solver,
                               localProblemFunctionToolkit->copy_local_reaction, contact, problem,
                               localproblem, reaction, local_opts, localreaction);

            if (options->iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION] ==
                SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE)
              localProblemFunctionToolkit->perform_relaxation(localreaction, &reaction[contact * 3], omega);

            light_error_2 = localProblemFunctionToolkit->light_error_squared(localreaction, &reaction[contact * 3]);
            light_error_sum += light_error_2;

            /* int test =100; */
            /* if (contact == test) */
            /* { */
            /*   printf("reaction[%i] = %16.8e\t",3*contact-1,reaction[3*contact]); */
            /*   printf("localreaction[%i] = %16.8e\n",2,localreaction[0]); */
            /* } */

            if (options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
              double squared_norm_localreaction = squared_norm(localreaction);
              int relative_convergence_criteria =
                  light_error_2 <= tmp_criteria1 * squared_norm_localreaction;
              int small_reaction_criteria = squared_norm_localreaction <= tmp_criteria2;

              if ((relative_convergence_criteria || small_reaction_criteria) && iter >= 10)
              /* if ((light_error_2 *squared_norm(localreaction) <=
               * tolerance*tolerance/(nc*nc*10)
               */
              /*      || squared_norm(localreaction) <=  (*norm_r* *norm_r/(nc*nc*1000))) */
              /*     && iter >=10) */
              {
                /* we  freeze the contact for n iterations*/
                freeze_contacts[contact] = options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT];

                DEBUG_EXPR(
                    printf("first criteria : light_error_2*squared_norm(localreaction) <= "
                           "tolerance*tolerance/(nc*nc*10) ==> %e <= %e, bool =%i\n",
                           light_error_2 * squared_norm(localreaction),
                           tolerance * tolerance / (nc * nc * 10),
                           relative_convergence_criteria);
                    printf("second criteria :  squared_norm(localreaction) <=  (*norm_r* "
                           "*norm_r/(nc*nc))/1000. ==> %e <= %e, bool =%i \n",
                           squared_norm(localreaction), *norm_r * *norm_r / (nc * nc * 1000),
                           small_reaction_criteria);
                    printf("Contact % i is freezed for %i steps\n", contact,
                           options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT]););
              }
            }

            if (options->iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] ==
                SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE)
              acceptLocalReactionFiltered(localproblem, local_opts, contact, iter,
                                          reaction, localreaction);
            else
              acceptLocalReactionUnconditionally(contact, reaction, localreaction);
          }
        }

        /* DEBUG_EXPR( */
        /*   if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0) */
        /*   { */
        /*     int frozen_contact=0; */
        /*     for(unsigned int ii = 0 ; ii < nc ; ++ii) if (freeze_contacts[ii] >0)
         * frozen_contact++; */
        /*     numerics_printf_verbose(1,"number of frozen contacts %i at iter : %i",
         * frozen_contact, iter ); */
        /*   } */
        /*   ); */

#pragma omp single
        {
          if (options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
              SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT) {
            error = calculateLightError(light_error_sum, nc, reaction, norm_r);
            hasNotConverged = determine_convergence(error, tolerance, iter, options);
          } else if (options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                     SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
            error = calculateLightError(light_error_sum, nc, reaction, norm_r);

            hasNotConverged = determine_convergence_with_full_final(
                problem, options, computeError, reaction, velocity, &tolerance, norm_q, error,
                iter);

            if (!(tolerance > 0.0)) {
              numerics_warning("fc3d_nsgs", "tolerance has to be positive!!");
              numerics_warning("fc3d_nsgs", "we stop the iterations");
              // I had to remove the break because OpenMP does not allow it.
              // break;
            }

          } else if (options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                     SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL) {
            error = calculateFullErrorAdaptiveInterval(problem, computeError, options, iter,
                                                       reaction, velocity, tolerance, norm_q);
            hasNotConverged = determine_convergence(error, tolerance, iter, options);
          }

          statsIterationCallback(problem, options, reaction, velocity, error);

          ++iter;
          light_error_sum = 0.0;
        }

        /* if(iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] >0) */
        /* { */
        /*   int frozen_contact=0; */
        /*   for(unsigned int i = 0 ; i < nc ; ++i) */
        /*   { */
        /*     if (freeze_contacts[i] >0) */
        /*     { */
        /*       frozen_contact++; */
        /*     } */
        /*   } */
        /*   printf("number of frozen contacts %i at iter : %i over number of contacts: %i\n",
         * frozen_contact, iter, nc ); */
        /* } */
      }

      // Clear memory
      solver_options_delete(local_opts);
      fc3d_local_problem_free(localproblem, problem);
      free(local_opts);
    }
  }

  /* Full criterium */
  if (options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
      SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
    error = calculateFullErrorFinal(problem, options, computeError, reaction, velocity,
                                    tolerance, norm_q);

    hasNotConverged = determine_convergence(error, SOLVER_TOL(options), iter, options);
  }

  *info = hasNotConverged;

  /** return parameter values */
  /* SOLVER_TOL(options) = tolerance; */
  SET_SOLVER_RESIDUAL(options, error);
  SET_SOLVER_ITER_DONE(options, iter);

  /* Restore original local solver tolerance */
  tolerance_manager_restore_local(&tol_manager, local_opts);

  /** Free memory **/

  if (problem->M->storageType == NM_SPARSE) {
    SBM_clear_block(problem->M->matrix1);
    SBM_clear(problem->M->matrix1);
    problem->M->matrix1 = matrix1;
  }
  localProblemFunctionToolkit->free_local_solver(problem, localproblem, local_opts);

  fc3d_local_problem_free(localproblem, problem);

  if (scontacts) free(scontacts);

  free(partition_size);
  for (size_t i = 0; i < n_colors; i++) free(partitions[i]);
  free(partitions);
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers FC3D_NSGS in the global solver registry, enabling:
 * - Dynamic solver lookup by ID
 * - Runtime solver introspection
 * - Elimination of giant switch statements in drivers
 */

static int fc3d_nsgs_graph_init_wrap(void* problem, SolverOptions* options) {
  /* set_default already called by solver_options_create */
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int fc3d_nsgs_graph_solve_wrap(void* problem, double* reaction,
                                double* velocity, SolverOptions* options) {
  int info = NUMERICS_OK;
  fc3d_nsgs((FrictionContactProblem*)problem, reaction, velocity, &info, options);
  return info;
}

static void fc3d_nsgs_graph_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_FRICTION_3D_NSGS_GRAPH, "SICONOS_FRICTION_3D_NSGS_GRAPH",
                "Non-smooth Gauss-Seidel for 3D Friction Contact",
                fc3d_nsgs_graph_init_wrap,
                fc3d_nsgs_graph_solve_wrap,
                fc3d_nsgs_graph_free_wrap,
                NULL,  /* error function */
                fc3d_nsgs_set_default,  /* set_default */
                1000,  /* default_max_iter */
                1e-4,  /* default_tol */
                0      /* is_local_solver */)