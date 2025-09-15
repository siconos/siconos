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

#include "FrictionContactProblem.h"  // for FrictionContac...
#include "Friction_cst.h"            // for SICONOS_FRICTI...
#include "NumericsArrays.h"          // for uint_shuffle
#include "NumericsFwd.h"             // for SolverOptions
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "SiconosBlas.h"    // for cblas_dnrm2
#include "SolverOptions.h"  // for SolverOptions
#include "SparseBlockMatrix.h"
#include "fc3d_2NCP_Glocker.h"                         // for NCPGlocker_update
#include "fc3d_NCPGlockerFixedPoint.h"                 // for fc3d_FixedP_in...
#include "fc3d_Path.h"                                 // for fc3d_Path_init...
#include "fc3d_Solvers.h"                              // for ComputeErrorPtr
#include "fc3d_compute_error.h"                        // for fc3d_compute_e...
#include "fc3d_local_problem_tools.h"                  // for fc3d_local_pro...
#include "fc3d_onecontact_nonsmooth_Newton_solvers.h"  // for fc3d_onecontac...
#include "fc3d_projection.h"                           // for fc3d_projectio...
#include "fc3d_unitary_enumerative.h"                  // for fc3d_unitary_e...
#include "numerics_verbose.h"                          // for numerics_printf
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */

#include "siconos_debug.h"                                     // for DEBUG_EXPR
#include "NumericsVector.h"

#include "gfc3d_ipm.h"
#include "projectionOnCone.h"        // for projectionOnCone
#include <time.h>


#include "op3x3.h"
#include "siconos_debug.h"  // for DEBUG_EXPR


//#define FCLIB_OUTPUT

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

static inline void performRelaxation_3(double localreaction[3], double *oldreaction,
                                       double omega) {
  localreaction[0] = omega * localreaction[0] + (1.0 - omega) * oldreaction[0];
  localreaction[1] = omega * localreaction[1] + (1.0 - omega) * oldreaction[1];
  localreaction[2] = omega * localreaction[2] + (1.0 - omega) * oldreaction[2];
}

static inline double light_error_squared_3(double localreaction[3], double *oldreaction) {
  return (pow(oldreaction[0] - localreaction[0], 2) +
          pow(oldreaction[1] - localreaction[1], 2) +
          pow(oldreaction[2] - localreaction[2], 2));
}

static inline double squared_norm_3(double localreaction[3]) {
  return (pow(localreaction[0], 2) + pow(localreaction[1], 2) + pow(localreaction[2], 2));
}

static double calculateLightError(double light_error_sum, unsigned int nc, double *reaction,
                                  double *norm_r) {
  double error = sqrt(light_error_sum);
  *norm_r = cblas_dnrm2(nc * 3, reaction, 1);
  if (fabs(*norm_r) > DBL_EPSILON) error /= (*norm_r);
  return error;
}

static void acceptLocalReactionUnconditionally(unsigned int contact, double *reaction,
                                               double localreaction[3]) {
  memcpy(&reaction[contact * 3], localreaction, sizeof(double) * 3);
}

static void statsIterationCallback(FrictionContactProblem *problem, SolverOptions *options,
                                   double *reaction, double *velocity, double error) {
  if (options->callback) {
    options->callback->collectStatsIteration(options->callback->env,
                                             problem->numberOfContacts * 3, reaction, velocity,
                                             error, NULL);
  }
}

static void fc3d_nsgs_update(int contact, FrictionContactProblem *problem,
                             FrictionContactProblem *localproblem, double *reaction,
                             SolverOptions *options) {
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
  fc3d_local_problem_compute_q(problem, localproblem, reaction, contact);

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];
}

void fc3d_nsgs_initialize_local_solver(
    struct LocalProblemFunctionToolkit *local_function_toolkit, ComputeErrorPtr *computeError,
    FrictionContactProblem *problem, FrictionContactProblem *localproblem,
    SolverOptions *options) {
  SolverOptions *localsolver_options = options->internalSolvers[0];

  *computeError = (ComputeErrorPtr)&fc3d_compute_error;

  if (problem->dimension == 3) {
    local_function_toolkit->copy_local_reaction = cpy3;
    local_function_toolkit->perform_relaxation = &performRelaxation_3;
    local_function_toolkit->light_error_squared = &light_error_squared_3;
    local_function_toolkit->squared_norm = &squared_norm_3;

  }

  /** Connect to local solver */
  switch (localsolver_options->solverId) {
    /* Projection */
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization: {
      local_function_toolkit->local_solver = &fc3d_projectionWithDiagonalization_solve;
      local_function_toolkit->update_local_problem =
          &fc3d_projectionWithDiagonalization_update;
      local_function_toolkit->free_local_solver = &fc3d_projection_free;
      fc3d_projection_initialize(problem, localproblem);
      break;
    }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone: {
      local_function_toolkit->local_solver = &fc3d_projectionOnCone_solve;
      local_function_toolkit->update_local_problem = &fc3d_nsgs_update;
      local_function_toolkit->free_local_solver = &fc3d_projection_free;
      fc3d_projection_initialize(problem, localproblem);
      break;
    }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration: {
      local_function_toolkit->local_solver = &fc3d_projectionOnConeWithLocalIteration_solve;
      local_function_toolkit->update_local_problem = &fc3d_nsgs_update;
      local_function_toolkit->free_local_solver =
          &fc3d_projectionOnConeWithLocalIteration_free;
      fc3d_projectionOnConeWithLocalIteration_initialize(problem, localproblem,
                                                         localsolver_options);
      break;
    }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization: {
      local_function_toolkit->local_solver = &fc3d_projectionOnCone_solve;
      local_function_toolkit->update_local_problem =
          &fc3d_projection_update_with_regularization;
      local_function_toolkit->free_local_solver = &fc3d_projection_with_regularization_free;
      fc3d_projection_initialize_with_regularization(problem, localproblem);
      break;
    }
    /* Newton solver (Alart-Curnier) */
    case SICONOS_FRICTION_3D_ONECONTACT_NSN: {
      local_function_toolkit->local_solver = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      local_function_toolkit->update_local_problem =
          &fc3d_onecontact_nonsmooth_Newton_AC_update;
      local_function_toolkit->free_local_solver =
          &fc3d_onecontact_nonsmooth_Newton_solvers_free;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem,
                                                          localsolver_options);
      break;
    }
    case SICONOS_FRICTION_3D_ONECONTACT_NSN_GP: {
      local_function_toolkit->local_solver = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      local_function_toolkit->update_local_problem =
          &fc3d_onecontact_nonsmooth_Newton_AC_update;
      local_function_toolkit->free_local_solver =
          &fc3d_onecontact_nonsmooth_Newton_solvers_free;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem,
                                                          localsolver_options);
      break;
    }
    case SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID: {
      local_function_toolkit->local_solver = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      local_function_toolkit->update_local_problem =
          &fc3d_onecontact_nonsmooth_Newton_AC_update;
      local_function_toolkit->free_local_solver =
          &fc3d_onecontact_nonsmooth_Newton_solvers_free;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem,
                                                          localsolver_options);
      break;
    } /* Newton solver (Glocker-Fischer-Burmeister)*/
    case SICONOS_FRICTION_3D_NCPGlockerFBNewton: {
      local_function_toolkit->local_solver = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      local_function_toolkit->update_local_problem = &NCPGlocker_update;
      local_function_toolkit->free_local_solver =
          &fc3d_onecontact_nonsmooth_Newton_solvers_free;
      // *computeError = &fake_compute_error;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem,
                                                          localsolver_options);
      break;
    }
    /* Path solver (Glocker Formulation) */
    case SICONOS_FRICTION_3D_NCPGlockerFBPATH: {
      local_function_toolkit->local_solver = &fc3d_Path_solve;
      local_function_toolkit->free_local_solver = (FreeLocalSolverPtr)&fc3d_Path_free;
      local_function_toolkit->update_local_problem = &NCPGlocker_update;

      // *computeError = &fake_compute_error;
      fc3d_Path_initialize(problem, localproblem, localsolver_options);
      break;
    }

    /* Fixed Point solver (Glocker Formulation) */
    case SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint: {
      local_function_toolkit->local_solver = &fc3d_FixedP_solve;
      local_function_toolkit->update_local_problem = &NCPGlocker_update;
      local_function_toolkit->free_local_solver = (FreeLocalSolverPtr)&fc3d_FixedP_free;
      /* *computeError = &fake_compute_error_nsgs; */
      fc3d_FixedP_initialize(problem, localproblem, localsolver_options);
      break;
    }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder: {
      local_function_toolkit->local_solver = &fc3d_projectionOnCylinder_solve;
      local_function_toolkit->update_local_problem = &fc3d_projectionOnCylinder_update;
      local_function_toolkit->free_local_solver =
          (FreeLocalSolverPtr)&fc3d_projectionOnCylinder_free;
      *computeError = (ComputeErrorPtr)&fc3d_Tresca_compute_error;
      fc3d_projectionOnCylinder_initialize(problem, localproblem, options);
      break;
    }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration: {
      local_function_toolkit->local_solver =
          &fc3d_projectionOnCylinderWithLocalIteration_solve;
      local_function_toolkit->update_local_problem = &fc3d_projectionOnCylinder_update;
      local_function_toolkit->free_local_solver =
          (FreeLocalSolverPtr)&fc3d_projectionOnCylinderWithLocalIteration_free;
      *computeError = (ComputeErrorPtr)&fc3d_Tresca_compute_error;
      fc3d_projectionOnCylinderWithLocalIteration_initialize(problem, localproblem, options,
                                                             localsolver_options);
      break;
    }
    case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC: {
      local_function_toolkit->local_solver = &fc3d_unitary_enumerative_solve;
      local_function_toolkit->update_local_problem = &fc3d_nsgs_update;
      local_function_toolkit->free_local_solver =
          (FreeLocalSolverPtr)&fc3d_unitary_enumerative_free;
      fc3d_unitary_enumerative_initialize(localproblem);
      break;
    }
    case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU: {
      local_function_toolkit->local_solver = &fc3d_unitary_enumerative_solve;
      local_function_toolkit->update_local_problem = &fc3d_nsgs_update;
      local_function_toolkit->free_local_solver =
          (FreeLocalSolverPtr)&fc3d_unitary_enumerative_free;
      fc3d_unitary_enumerative_initialize(localproblem);
      break;
    }
  default: {
      numerics_error("fc3d_nsgs_initialize_local_solver",
                     "Numerics, fc3d_nsgs failed. Unknown internal solver : %s.\n",
                     solver_options_id_to_name(localsolver_options->solverId));
    }
  }
}

static unsigned int *allocShuffledContacts(FrictionContactProblem *problem,
                                           SolverOptions *options) {
  unsigned int *scontacts = 0;
  unsigned int nc = problem->numberOfContacts;
  if (options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
          SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE ||
      options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
          SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP) {
    if (options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED] > 0) {
      srand((unsigned int)options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED]);
    } else
      srand(1);
    scontacts = (unsigned int *)malloc(nc * sizeof(unsigned int));
    for (unsigned int i = 0; i < nc; ++i) {
      scontacts[i] = i;
    }
    uint_shuffle(scontacts, nc);
  }
  return scontacts;
}
static unsigned int *allocfreezingContacts(FrictionContactProblem *problem,
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

static int solveLocalReaction(UpdatePtr update_localproblem, SolverPtr local_solver,
                              CopyLocalReactionPtr copyLocalReaction, unsigned int contact,
                              FrictionContactProblem *problem,
                              FrictionContactProblem *localproblem, double *reaction,
                              SolverOptions *localsolver_options, double localreaction[3]) {
  (*update_localproblem)(contact, problem, localproblem, reaction, localsolver_options);

  localsolver_options->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER] = contact;

  copyLocalReaction(&(reaction[contact * problem->dimension]), localreaction);

  return (*local_solver)(localproblem, localreaction, localsolver_options);
}

int file_exists(const char *fname) {
  FILE *file;
  if ((file = fopen(fname, "r"))) {
    fclose(file);
    return 1;
  }
  return 0;
}

static void acceptLocalReactionFiltered(FrictionContactProblem *localproblem,
                                        SolverOptions *localsolver_options,
                                        unsigned int contact, unsigned int iter,
                                        double *reaction, double localreaction[3]) {
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
      char *title = (char *)malloc(n * sizeof(char));
      strcpy(title, "Bad local problem dump in hdf5");
      char *description = (char *)malloc(n * sizeof(char));
      strcpy(description, "Rewriting in hdf5 from siconos ");
      strcat(description, fname);
      strcat(description, " in FCLIB format");
      char *mathInfo = (char *)malloc(n * sizeof(char));
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
    memcpy(&reaction[contact * localproblem->dimension], localreaction,
           sizeof(double) * localproblem->dimension);

}

static double calculateFullErrorAdaptiveInterval(FrictionContactProblem *problem,
                                                 ComputeErrorPtr computeError,
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
        "--------------- FC3D - NSGS - Iteration %i "
        "options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = %i, "
        "options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] = % i",
        iter, options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY],
        options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION]);
  } else
    (*computeError)(problem, reaction, velocity, tolerance, options, norm_q, &error);

  return error;
}

static double calculateFullErrorFinal(FrictionContactProblem *problem, SolverOptions *options,
                                      ComputeErrorPtr computeError, double *reaction,
                                      double *velocity, double tolerance, double norm_q) {
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
                                 SolverOptions *options) {
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

static int determine_convergence_with_full_final(FrictionContactProblem *problem,
                                                 SolverOptions *options,
                                                 ComputeErrorPtr computeError,
                                                 double *reaction, double *velocity,
                                                 double *tolerance, double norm_q,
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

void fc3d_nsgs(FrictionContactProblem *problem, double *reaction, double *velocity, int *info,
               SolverOptions *options) {
  /* verbose=1; */

  /* int and double parameters */
  int *iparam = options->iparam;
  double *dparam = options->dparam;

  /* Number of contacts */
  unsigned int nc = problem->numberOfContacts;

  /* Maximum number of iters */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];

  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  double norm_q = cblas_dnrm2(nc * 3, problem->q, 1);
  double omega = dparam[SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE];

  double norm_r[] = {1e24};
  if (options->numberOfInternalSolvers < 1) {
    numerics_error("fc3d_nsgs",
                   "The NSGS method needs options for the internal solvers, "
                   "options[0].numberOfInternalSolvers should be >= 1");
  }

  SolverOptions *localsolver_options = options->internalSolvers[0];
  ComputeErrorPtr computeError = NULL;

  struct LocalProblemFunctionToolkit *localProblemFunctionToolkit =
      localProblemFunctionToolkit_new();
  /* localProblemFunctionToolkit_display(localProblemFunctionToolkit); */

  FrictionContactProblem *localproblem;

  double localreaction[3];

  /*****  NSGS Iterations *****/
  int iter = 0;      /* Current iteration number */

  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  unsigned int contact; /* Number of the current row of blocks in M */
  unsigned int *scontacts = NULL;
  unsigned int *freeze_contacts = NULL;

  if (*info == 0) return;

  SparseBlockStructuredMatrix *matrix1 = problem->M->matrix1;
  if (problem->M->storageType == NM_SPARSE) {
    if (problem->M->matrix1) {
      printf("Warning matrix 1 different from NULL");
    }

    problem->M->matrix1 = NM_extract_diagonal_blocks(problem->M, problem->dimension);
  }

  /*****  Initialize various solver options *****/
  localproblem = fc3d_local_problem_allocate(problem);

  fc3d_nsgs_initialize_local_solver(localProblemFunctionToolkit, &computeError, problem,
                                    localproblem, options);

  /* localProblemFunctionToolkit_display(localProblemFunctionToolkit); */
  scontacts = allocShuffledContacts(problem, options);
  freeze_contacts = allocfreezingContacts(problem, options);
  /*****  Check solver options *****/
  if (!(iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE ||
        iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE ||
        iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
            SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP)) {
    numerics_error("fc3d_nsgs",
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
    numerics_error("fc3d_nsgs",
                   "iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] must be equal to "
                   "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL (0), "
                   "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL (1), "
                   "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT (2) or "
                   "SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE (3)");
    return;
  }
// FILE *iterates = NULL;
  /*****  NSGS Iterations *****/

  /* A special case for the most common options (should correspond
   * with mechanics_run.py **/
  if (iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] == SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE &&
      iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] == 0 &&
      iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION] ==
          SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE &&
      iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] ==
          SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE &&
      iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
          SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT) {
    while ((iter < itermax) && (hasNotConverged > 0)) {
      ++iter;
      double light_error_sum = 0.0;

      fc3d_set_internalsolver_tolerance(problem, options, localsolver_options, error);

      for (unsigned int i = 0; i < nc; ++i) {
        contact = i;

        solveLocalReaction(localProblemFunctionToolkit->update_local_problem,
                           localProblemFunctionToolkit->local_solver,
                           localProblemFunctionToolkit->copy_local_reaction, contact, problem,
                           localproblem, reaction, localsolver_options, localreaction);

        light_error_sum += localProblemFunctionToolkit->light_error_squared(
            localreaction, &reaction[contact * 3]);

        /* #if 0 */
        acceptLocalReactionFiltered(localproblem, localsolver_options, contact, iter, reaction,
                                    localreaction);
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
    /* verbose=1; */
    while ((iter < itermax) && (hasNotConverged > 0)) {
      ++iter;
      double light_error_sum = 0.0;
      double light_error_2 = 0.0;
      fc3d_set_internalsolver_tolerance(problem, options, localsolver_options, error);

      unsigned int number_of_freezed_contact = 0;
      double tmp_criteria1 = tolerance * tolerance * 100 * 100;
      double tmp_criteria2 = *norm_r * *norm_r / (nc * nc * 1000);

      if (iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
        for (unsigned int i = 0; i < nc; ++i) {
          if (freeze_contacts[i] > 0) number_of_freezed_contact++;
        }
        if (number_of_freezed_contact >= nc - 1) {
          // printf("number of freezed contact too large\n");
          for (unsigned int c = 0; c < nc; ++c) freeze_contacts[c] = 0;
        }
      }
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
            continue;
          }
        }

        solveLocalReaction(localProblemFunctionToolkit->update_local_problem,
                           localProblemFunctionToolkit->local_solver,
                           localProblemFunctionToolkit->copy_local_reaction, contact, problem,
                           localproblem, reaction, localsolver_options, localreaction);

        if (iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION] ==
            SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE)
          localProblemFunctionToolkit->perform_relaxation(localreaction,
                                                          &reaction[contact * 3], omega);

        light_error_2 = localProblemFunctionToolkit->light_error_squared(
            localreaction, &reaction[contact * 3]);

        light_error_sum += light_error_2;

        /* int test =100; */
        /* if (contact == test) */
        /* { */
        /*   printf("reaction[%i] = %16.8e\t",3*contact-1,reaction[3*contact]); */
        /*   printf("localreaction[%i] = %16.8e\n",2,localreaction[0]); */
        /* } */

        if (iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
          double squared_norm_localreaction =
              localProblemFunctionToolkit->squared_norm(localreaction);
          int relative_convergence_criteria =
              light_error_2 <= tmp_criteria1 * squared_norm_localreaction;
          int small_reaction_criteria = squared_norm_localreaction <= tmp_criteria2;

          if ((relative_convergence_criteria || small_reaction_criteria) && iter >= 10)
          /* if ((light_error_2 *squared_norm(localreaction) <= tolerance*tolerance/(nc*nc*10)
           */
          /*      || squared_norm(localreaction) <=  (*norm_r* *norm_r/(nc*nc*1000))) */
          /*     && iter >=10) */
          {
            /* we  freeze the contact for n iterations*/

            freeze_contacts[contact] = iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT];

            DEBUG_EXPR(printf("first criteria : light_error_2*squared_norm(localreaction) <= "
                              "tolerance*tolerance/(nc*nc*10) ==> %e <= %e, bool =%i\n",
                              light_error_2 * squared_norm(localreaction),
                              tolerance * tolerance / (nc * nc * 10),
                              relative_convergence_criteria);
                       printf("second criteria :  squared_norm(localreaction) <=  (*norm_r* "
                              "*norm_r/(nc*nc))/1000. ==> %e <= %e, bool =%i \n",
                              squared_norm(localreaction),
                              *norm_r * *norm_r / (nc * nc * 1000), small_reaction_criteria);
                       printf("Contact % i is freezed for %i steps\n", contact,
                              iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT]););
          }
        }

        if (iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] ==
            SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE)
          acceptLocalReactionFiltered(localproblem, localsolver_options, contact, iter,
                                      reaction, localreaction);
        else
          acceptLocalReactionUnconditionally(contact, reaction, localreaction);
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

        if (!(tolerance > 0.0)) {
          numerics_warning("fc3d_nsgs", "tolerance has to be positive!!");
          numerics_warning("fc3d_nsgs", "we stop the iters");
          break;
        }

      } else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                 SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL) {
        error = calculateFullErrorAdaptiveInterval(problem, computeError, options, iter,
                                                   reaction, velocity, tolerance, norm_q);
        hasNotConverged = determine_convergence(error, tolerance, iter, options);
      }

      statsIterationCallback(problem, options, reaction, velocity, error);

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
  }
  /* Full criterium */
  if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
      SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
    error = calculateFullErrorFinal(problem, options, computeError, reaction, velocity,
                                    tolerance, norm_q);

    hasNotConverged = determine_convergence(error, dparam[SICONOS_DPARAM_TOL], iter, options);
  }


  *info = hasNotConverged;

  /** return parameter values */
  /* dparam[SICONOS_DPARAM_TOL] = tolerance; */
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

  /** Free memory **/

  if (problem->M->storageType == NM_SPARSE) {
    SBM_clear_block(problem->M->matrix1);
    SBM_clear(problem->M->matrix1);
    problem->M->matrix1 = matrix1;
  }
  localProblemFunctionToolkit->free_local_solver(problem, localproblem, localsolver_options);

  fc3d_local_problem_free(localproblem, problem);

  if (scontacts) free(scontacts);
   
}


void fc3d_nsgs_set_default(SolverOptions *options) {
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
  // Internal solver
  assert(options->numberOfInternalSolvers == 1);
  options->internalSolvers[0] =
      solver_options_create(SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID);
  // Printing in the same style as in IPM solver
  options->iparam[SICONOS_FRICTION_3D_NSGS_PRINTING_LIKE_IPM] = SICONOS_FRICTION_3D_NSGS_PRINTING_LIKE_IPM_TRUE;

}
