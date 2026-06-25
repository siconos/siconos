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
#include "cohesive_friction_3d_nsgs.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "CohesiveFrictionContactProblem.h"
#include "CohesiveFrictionContact_options.h"
#include "FrictionContactProblem.h"
#include "NonSmoothGaussSeidel_options.h"
#include "NumericsMatrix.h"
#include "NumericsVector.h"
#include "SiconosBlas.h"
#include "SolverOptions.h"
#include "cohesive_friction_3d_compute_error.h"
#include "cohesive_friction_3d_local_problem_tools.h"
#include "cohesive_friction_3d_projection.h"
#include "fc3d_projection.h"
#include "fc3d_onecontact_nonsmooth_Newton_solvers.h"
#include "fc3d_short_names.h"
#include "naming_conventions.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "NumericsArrays.h"
#include "op3x3.h"
#include "solver_registry.h"
#include "tolerance_manager.h"
#include "Friction_tools.h"
#include "nsgs_generic.h"

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"

/** pointer to function used to update velocity and compute error */
typedef int (*CohesiveComputeErrorPtr)(CohesiveFrictionContactProblem *, double *, double *, double,
				       SolverOptions *, double, double *);


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

static double calculateLightError(double light_error_sum, unsigned int nc, double* reaction,
                                  double* norm_r) {
  double error = sqrt(light_error_sum);
  *norm_r = cblas_dnrm2(nc * 3, reaction, 1);
  if (fabs(*norm_r) > DBL_EPSILON) error /= (*norm_r);
  return error;
}

static void acceptLocalReactionUnconditionally(unsigned int contact, double* reaction,
                                               double localreaction[3]) {
  memcpy(&reaction[contact * 3], localreaction, sizeof(double) * 3);
}

static void statsIterationCallback(CohesiveFrictionContactProblem* problem,
                                   SolverOptions* options, double* reaction, double* velocity,
                                   double error) {
  if (options->callback) {
    options->callback->collectStatsIteration(options->callback->env,
                                             problem->numberOfContacts * 3, reaction, velocity,
                                             error, NULL);
  }
}

static void cohesive_fc3d_nsgs_update(int contact, CohesiveFrictionContactProblem* problem,
                                      CohesiveFrictionContactProblem* localproblem, double* reaction,
                                      SolverOptions* options) {
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */
  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  cohesive_friction_3d_local_problem_fill_M(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products
     MLocal.reactionBlock, excluding the block corresponding to the current contact. ****/
  cohesive_friction_3d_local_problem_compute_q(problem, localproblem, reaction, contact);

  /* Friction coefficient for current block*/
  
  int nc = problem->numberOfContacts;
  if (contact < nc) {
    localproblem->mu[0] = problem->mu[contact];
    localproblem->c_n[0] = 0.0;
    localproblem->c_t[0] = 0.0;
  } else {
    localproblem->mu[0] = 0.0;
    localproblem->c_n[0] = problem->c_n[contact-nc];
    localproblem->c_t[0] = problem->c_t[contact-nc];

    }  
}

static int cohesive_fc3d_nsgs_initialize_local_solver(
    struct CohesiveLocalProblemFunctionToolkit* local_function_toolkit,
    CohesiveFrictionContactProblem* problem,
    FrictionContactProblem* localproblem_contact,
    CohesiveFrictionContactProblem* localproblem_cohesion,    
    SolverOptions* options) {

  

  SolverOptions* local_opts_contact = options->internalSolvers[0];
  SolverOptions* local_opts_cohesion = options->internalSolvers[1];

  if (problem->dimension == 3) {
    local_function_toolkit->copy_local_reaction = cpy3;
    local_function_toolkit->perform_relaxation = &performRelaxation_3;
    local_function_toolkit->light_error_squared = &light_error_squared_3;
    local_function_toolkit->squared_norm = &squared_norm_3;
  }


  /** Create a Frictioncontactproblem for the initialization of contact local problem*/
  FrictionContactProblem* fc3d_problem = frictionContactProblem_new();
  fc3d_problem->numberOfContacts = problem->numberOfContacts;
  fc3d_problem->M = problem->M;

  
  /** Connect to local solver */
  switch (local_opts_contact->solverId) {
    /* Projection */
    case OC_PROJ: {
      local_function_toolkit->local_solver_contact = &fc3d_projectionOnCone_solve;
      local_function_toolkit->update_local_problem = &cohesive_fc3d_nsgs_update;
      local_function_toolkit->free_local_solver_contact = &fc3d_projection_free;
      fc3d_projection_initialize(fc3d_problem);
      break;
    }
    case OC_PROJ_LI: {
      local_function_toolkit->local_solver_contact = &fc3d_projectionOnConeWithLocalIteration_solve;
      local_function_toolkit->update_local_problem = &cohesive_fc3d_nsgs_update;
      local_function_toolkit->free_local_solver_contact =
          &fc3d_projectionOnConeWithLocalIteration_free;
      fc3d_projectionOnConeWithLocalIteration_initialize(fc3d_problem, local_opts_contact);
      break;
    }

    /* Newton solver (Alart-Curnier) */
    case OC_NSN: {
      local_function_toolkit->local_solver_contact = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      local_function_toolkit->update_local_problem =
	&cohesive_fc3d_nsgs_update;
      local_function_toolkit->free_local_solver_contact =
          &fc3d_onecontact_nonsmooth_Newton_solvers_free;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(fc3d_problem, local_opts_contact);
      break;
    }
    case OC_NSN_GP: {
      local_function_toolkit->local_solver_contact = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      local_function_toolkit->update_local_problem =
          &cohesive_fc3d_nsgs_update;
      local_function_toolkit->free_local_solver_contact =
          &fc3d_onecontact_nonsmooth_Newton_solvers_free;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(fc3d_problem, local_opts_contact);
      break;
    }
    case OC_NSN_GP_HYBRID: {
      local_function_toolkit->local_solver_contact = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      local_function_toolkit->update_local_problem =
          &cohesive_fc3d_nsgs_update;
      local_function_toolkit->free_local_solver_contact =
          &fc3d_onecontact_nonsmooth_Newton_solvers_free;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(fc3d_problem, local_opts_contact);
      break;
    }

    default: {
      return numerics_error(
          "cohesive_fc3d_nsgs_initialize_local_solver_contact",
          "Numerics, cohesive_fc3d_nsgs failed. Unknown internal solver : %s.\n",
          solver_options_id_to_name(local_opts_contact->solverId));
    }
  }
    /** Connect to local solver */
  switch (local_opts_cohesion->solverId) {
    /* Projection */
    case SICONOS_COHESIVE_FRICTION_3D_PROJECTION: {
      local_function_toolkit->local_solver_cohesion = &cohesive_friction_3d_projection_solve;
      local_function_toolkit->update_local_problem = &cohesive_fc3d_nsgs_update;
       local_function_toolkit->free_local_solver_cohesion = &cohesive_friction_3d_projection_free;
      cohesive_friction_3d_projection_initialize(problem);
      break;
    }
  default: {
      return numerics_error(
          "cohesive_fc3d_nsgs_initialize_local_solver",
          "Numerics, cohesive_fc3d_nsgs failed. Unknown internal solver : %s.\n",
          solver_options_id_to_name(local_opts_contact->solverId));
    }
  }
 
  
  return 0;
}
static unsigned int* allocShuffledContacts(CohesiveFrictionContactProblem* problem,
                                           SolverOptions* options) {
  unsigned int* scontacts = 0;
  unsigned int nc = problem->numberOfContacts;
  unsigned int ncoh = problem->numberOfCohesivePoints;  
  if (options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
          SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE ||
      options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] ==
          SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP) {
    if (options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED] > 0) {
      srand((unsigned int)options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED]);
    } else
      srand(1);
    scontacts = (unsigned int*)malloc((nc+ncoh) * sizeof(unsigned int));
    for (unsigned int i = 0; i < (nc+ncoh) ; ++i) {
      scontacts[i] = i;
    }
    uint_shuffle(scontacts, (nc+ncoh) );
  }
  return scontacts;
}
static unsigned int* allocfreezingContacts(CohesiveFrictionContactProblem* problem,
                                           SolverOptions* options) {
  unsigned int* fcontacts = 0;
  unsigned int nc = problem->numberOfContacts;
  unsigned int ncoh = problem->numberOfCohesivePoints;
  if (options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
    fcontacts = (unsigned int*)malloc((nc+ncoh) * sizeof(unsigned int));
    for (unsigned int i = 0; i < (nc+ncoh) ; ++i) {
      fcontacts[i] = 0;
    }
  }
  return fcontacts;
}

static int solveLocalReaction(
    struct CohesiveLocalProblemFunctionToolkit* localProblemFunctionToolkit,
    unsigned int contact, CohesiveFrictionContactProblem* problem,
    CohesiveFrictionContactProblem* localproblem, FrictionContactProblem* localproblem_contact,
    double* reaction, SolverOptions* local_opts, double localreaction[3]) {

  
  (*localProblemFunctionToolkit->update_local_problem)(contact, problem, localproblem, reaction, local_opts);

  local_opts->iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER] = contact;

  localProblemFunctionToolkit->copy_local_reaction(&(reaction[contact * problem->dimension]), localreaction);
  if (contact < problem->numberOfContacts) {

    localproblem_contact->M = localproblem->M;
    localproblem_contact->q = localproblem->q;
    localproblem_contact->mu = localproblem->mu;
    
    return (*localProblemFunctionToolkit->local_solver_contact)(localproblem_contact, localreaction,
                                                                local_opts);
    }    
  else {
    return (*localProblemFunctionToolkit->local_solver_cohesion)(localproblem, localreaction, local_opts);
  }
  //  return -1;  
}

static int file_exists(const char* fname) {
  FILE* file;
  if ((file = fopen(fname, "r"))) {
    fclose(file);
    return 1;
  }
  return 0;
}  
static void acceptLocalReactionFiltered(int dimension,
                                        SolverOptions* local_opts, unsigned int contact,
                                        unsigned int iter, double* reaction,
                                        double localreaction[3]) {
  if (isnan(SOLVER_RESIDUAL(local_opts)) || isinf(SOLVER_RESIDUAL(local_opts)) ||
      SOLVER_RESIDUAL(local_opts) > 1.0) {
    DEBUG_EXPR(frictionContact_display(localproblem));

    DEBUG_PRINTF(
        "Discard local reaction for contact %i at iteration %i "
        "with local_error = %e\n",
        contact, iter, SOLVER_RESIDUAL(local_opts));


    numerics_printf(
        "Discard local reaction for contact %i at iteration %i "
        "with local_error = %e",
        contact, iter, SOLVER_RESIDUAL(local_opts));
  } else
    memcpy(&reaction[contact * dimension], localreaction,
           sizeof(double) * dimension);
}

static double calculateFullErrorFinal(CohesiveFrictionContactProblem* problem, SolverOptions* options,
                                      CohesiveComputeErrorPtr computeError, double* reaction,
                                      double* velocity, double tolerance, double norm_q) {
  double absolute_error;
  (*computeError)(problem, reaction, velocity, tolerance, options, norm_q, &absolute_error);

  if (verbose > 0) {
    if (absolute_error > SOLVER_TOL(options)) {
      numerics_printf(
          "------- FC3D - NSGS - Warning absolute "
          "Residual = %14.7e is larger than required precision = %14.7e",
          absolute_error, SOLVER_TOL(options));
    } else {
      numerics_printf(
          "------- FC3D - NSGS - absolute "
          "Residual = %14.7e is smaller than required precision = %14.7e",
          absolute_error, SOLVER_TOL(options));
    }
  }
  return absolute_error;
}
/** Check convergence with full error verification and tolerance adaptation
 *
 * This function checks if the NSGS solver has converged by:
 * 1. Checking if incremental error is below working tolerance
 * 2. If yes, computing full error and checking against user tolerance
 * 3. Adapting tolerance if incremental converged but full didn't
 *
 * \param[in] problem Friction contact problem
 * \param[in] options Solver options
 * \param[in] computeError Function to compute full error
 * \param[in] reaction Current reaction vector
 * \param[in] velocity Current velocity vector
 * \param[in,out] tm Tolerance manager (handles adaptation)
 * \param[in] norm_q Norm of q vector
 * \param[in] incr_error Incremental error
 * \param[in] iter Current iteration
 * \return 0 if converged, 1 if not converged
 */
static int check_convergence_with_adaptation(CohesiveFrictionContactProblem* problem,
                                             SolverOptions* options,
                                             CohesiveComputeErrorPtr computeError, double* reaction,
                                             double* velocity, ToleranceManager* tm,
                                             double norm_q, double incr_error,
                                             double* full_error, int iter) {
  /* Check if incremental error is below working tolerance */
  if (incr_error >= tm->working_tolerance) {
    /* numerics_printf( */
    /*     "--------------- FC3D - NSGS - Iteration %i " */
    /*     "Residual = %14.7e > %7.3e", */
    /*     iter, incr_error, tm->working_tolerance); */
    return 1; /* Not converged */
  }

  /* Incremental error converged - check full error */
  /* numerics_printf( */
  /*     "--------------- FC3D - NSGS - Iteration %i " */
  /*     "Residual = %14.7e < %7.3e", */
  /*     iter, incr_error, tm->working_tolerance); */

  *full_error = calculateFullErrorFinal(problem, options, computeError, reaction, velocity,
                                        SOLVER_TOL(options), norm_q);

  /* Use tolerance manager to handle adaptation logic */
  SolverOptions* local_opts =
      (options->numberOfInternalSolvers > 0) ? options->internalSolvers[0] : NULL;
  int converged =
      tolerance_manager_check_convergence(tm, local_opts, *full_error, incr_error, verbose);

  if (converged == 0) {
    numerics_printf(
        "------- FC3D - NSGS - The incremental precision is sufficient to reach accuracy "
        "to %e",
        tm->working_tolerance);
  }

  return converged;
}

/* Deprecated: Use check_convergence_with_adaptation() with ToleranceManager */
static int determine_convergence_with_full_final(CohesiveFrictionContactProblem* problem,
                                                 SolverOptions* options,
                                                 CohesiveComputeErrorPtr computeError,
                                                 double* reaction, double* velocity,
                                                 double* tolerance, double norm_q,
                                                 double error, double* full_error, int iter) {
  /* Create temporary tolerance manager for backward compatibility */
  ToleranceManager tm;
  SolverOptions* local_opts =
      (options->numberOfInternalSolvers > 0) ? options->internalSolvers[0] : NULL;
  tolerance_manager_init(&tm, SOLVER_TOL(options), local_opts);
  tm.working_tolerance = *tolerance;

  int result =
      check_convergence_with_adaptation(problem, options, computeError, reaction, velocity,
                                        &tm, norm_q, error, full_error, iter);

  /* Sync back the adapted tolerance */
  *tolerance = tm.working_tolerance;

  return result;
}
static double calculateFullErrorAdaptiveInterval(CohesiveFrictionContactProblem* problem,
                                                 CohesiveComputeErrorPtr computeError,
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

/* ===========================================================================
 * Main NSGS Solver
 * =========================================================================== */

int cohesive_friction_3d_nsgs(CohesiveFrictionContactProblem* problem, double* reaction,
                              double* velocity, SolverOptions* options) {
  DEBUG_BEGIN("cohesive_friction_3d_nsgs(...)\n");

  if (!problem || !reaction || !velocity || !options) {
    return numerics_error("cohesive_friction_3d_nsgs", "NULL pointer argument");
  }
  /* Number of contacts */
  unsigned int nc = problem->numberOfContacts;
  unsigned int ncoh = problem->numberOfCohesivePoints;

  /* Maximum number of iterations */
  int itermax = SOLVER_MAX_ITER(options);

  /* Tolerance setup with unified tolerance manager */
  double norm_q = cblas_dnrm2((nc+ncoh) * 3, problem->q, 1);
  double omega = options->dparam[SICONOS_NSGS_RELAXATION_VALUE];

  double norm_r[] = {1e24};
  if (options->numberOfInternalSolvers < 2) {
    return numerics_error("cohesive_friction_3d_nsgs",
                          "The NSGS method needs options for the internal solvers, "
                          "options[0].numberOfInternalSolvers should be >= 2");
  }

  /* Get local solver options - use consistent naming */
  SolverOptions* local_opts_contact = options->internalSolvers[0];
  SolverOptions* local_opts_cohesion = options->internalSolvers[1];

  /* Initialize tolerance manager for unified tolerance handling */
  ToleranceManager tol_manager;
  tolerance_manager_init(&tol_manager, SOLVER_TOL(options), local_opts_contact);

  /* Working tolerance (may be adapted during iterations) */
  double tolerance = tol_manager.working_tolerance;

  CohesiveComputeErrorPtr computeError = NULL;
  computeError = &cohesive_friction_3d_compute_error;

  struct CohesiveLocalProblemFunctionToolkit* localProblemFunctionToolkit =
      cohesiveLocalProblemFunctionToolkit_new();
  
  /* localProblemFunctionToolkit_display(localProblemFunctionToolkit); */

  FrictionContactProblem* localproblem_contact;
  
  CohesiveFrictionContactProblem* localproblem_cohesion;
  
  double localreaction[3];

  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */

  double incr_error = INFINITY; /* Current error */
  double full_error = INFINITY; /* Current error */
  int hasNotConverged = 1;
  unsigned int contact; /* Number of the current row of blocks in M */
  unsigned int* scontacts = NULL;
  unsigned int* freeze_contacts = NULL;
  int frozen_contact = 0;

  SparseBlockStructuredMatrix* matrix1 = problem->M->matrix1;
  if (problem->M->storageType == NM_SPARSE) {
    if (problem->M->matrix1) {
      printf("Warning matrix 1 different from NULL");
    }

    problem->M->matrix1 = NM_extract_diagonal_blocks(problem->M, problem->dimension);
  }

  /*****  Initialize various solver options *****/
 
  localproblem_cohesion = cohesive_friction_3d_local_problem_allocate(problem->M->storageType);
  localproblem_contact = frictionContactProblem_new(); // wrap onto the local_problem_cohesion to call friction contact solver
  

  cohesive_fc3d_nsgs_initialize_local_solver(localProblemFunctionToolkit, problem,
                                             localproblem_contact, localproblem_cohesion, options);

  /* localProblemFunctionToolkit_display(localProblemFunctionToolkit); */
  scontacts = allocShuffledContacts(problem, options);
  freeze_contacts = allocfreezingContacts(problem, options);

  /*****  Check solver options *****/
  if (!(options->iparam[SICONOS_NSGS_SHUFFLE] == SICONOS_NSGS_SHUFFLE_FALSE ||
        options->iparam[SICONOS_NSGS_SHUFFLE] == SICONOS_NSGS_SHUFFLE_TRUE ||
        options->iparam[SICONOS_NSGS_SHUFFLE] == SICONOS_NSGS_SHUFFLE_EACH_LOOP)) {
    return  numerics_error("cohesive_fc3d_nsgs",
                           "options->iparam[SICONOS_NSGS_SHUFFLE] must be equal to "
                           "SICONOS_NSGS_SHUFFLE_FALSE (0), "
                           "SICONOS_NSGS_SHUFFLE_TRUE (1) or "
                           "SICONOS_NSGS_SHUFFLE_TRUE_EACH_LOOP (2)");
  }

  if (!(options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
            SICONOS_NSGS_ERROR_EVALUATION_FULL ||
        options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
            SICONOS_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL ||
        options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
            SICONOS_NSGS_ERROR_EVALUATION_LIGHT ||
        options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
            SICONOS_NSGS_ERROR_EVALUATION_ADAPTIVE)) {
    return numerics_error(
        "cohesive_fc3d_nsgs",
        "options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] must be equal to "
        "SICONOS_NSGS_ERROR_EVALUATION_FULL (0), "
        "SICONOS_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL (1), "
        "SICONOS_NSGS_ERROR_EVALUATION_LIGHT (2) or "
        "SICONOS_NSGS_ERROR_EVALUATION_ADAPTIVE (3)");

  }
  // FILE *iterates = NULL;
  /*****  NSGS Iterations *****/

  double* light_error_2 = light_error_2 = calloc(nc+ncoh, sizeof(double));

  //verbose=1;
  while ((iter < itermax) && (hasNotConverged > 0)) {
    ++iter;
    double light_error_sum = 0.0;

    fc3d_set_internalsolver_tolerance(nc, options, local_opts_contact, incr_error);

    unsigned int number_of_freezed_contact = 0;
    double tmp_criteria1 = tolerance * tolerance / (nc * nc * 1000);
    double tmp_criteria2 = *norm_r * *norm_r / (nc * nc * 1000);

    if (options->iparam[SICONOS_NSGS_FREEZING_CONTACT] > 0) {
      for (unsigned int i = 0; i < nc; ++i) {
        if (freeze_contacts[i] > 0) number_of_freezed_contact++;
      }

      if (number_of_freezed_contact >= nc - 1) {
        numerics_printf_verbose(1, "number of freezed contact is too large : %i\n",
                                number_of_freezed_contact);
        for (unsigned int c = 0; c < nc; ++c) freeze_contacts[c] = 0;
      }
    }
    for (unsigned int i = 0; i < nc+ncoh; ++i) {
      if (options->iparam[SICONOS_NSGS_SHUFFLE] == SICONOS_NSGS_SHUFFLE_TRUE ||
          options->iparam[SICONOS_NSGS_SHUFFLE] == SICONOS_NSGS_SHUFFLE_EACH_LOOP) {
        if (options->iparam[SICONOS_NSGS_SHUFFLE] == SICONOS_NSGS_SHUFFLE_EACH_LOOP)
          uint_shuffle(scontacts, nc);
        contact = scontacts[i];
      } else
        contact = i;

      /* if (options->iparam[SICONOS_NSGS_FREEZING_CONTACT] > 0) { */
      /*   if (freeze_contacts[contact] > 0) { */
      /*     /\* we skip freeze contacts *\/ */
      /*     freeze_contacts[contact] -= 1; */
      /*     light_error_sum += light_error_2[contact]; */
      /*     continue; */
      /*   } */
      /* } */

      solveLocalReaction(localProblemFunctionToolkit, contact, problem,
                         localproblem_cohesion, localproblem_contact, reaction, local_opts_contact, localreaction);

      if (options->iparam[SICONOS_NSGS_RELAXATION] == SICONOS_NSGS_RELAXATION_TRUE)
        localProblemFunctionToolkit->perform_relaxation(localreaction, &reaction[contact * 3],
                                                        omega);

      light_error_2[contact] = localProblemFunctionToolkit->light_error_squared(
          localreaction, &reaction[contact * 3]);

      light_error_sum += light_error_2[contact];

      /* int test =100; */
      /* if (contact == test) */
      /* { */
      /*   printf("reaction[%i] = %16.8e\t",3*contact-1,reaction[3*contact]); */
      /*   printf("localreaction[%i] = %16.8e\n",2,localreaction[0]); */
      /* } */

      /* if (options->iparam[SICONOS_NSGS_FREEZING_CONTACT] > 0) { */
      /*   double squared_norm_localreaction = */
      /*       localProblemFunctionToolkit->squared_norm(localreaction); */
      /*   int relative_convergence_criteria = */
      /*       light_error_2[contact] <= tmp_criteria1 * squared_norm_localreaction; */
      /*   int small_reaction_criteria = squared_norm_localreaction <= tmp_criteria2; */

      /*   if ((relative_convergence_criteria || small_reaction_criteria) && iter >= 10) */
      /*   /\* if ((light_error_2 *squared_norm(localreaction) <= tolerance*tolerance/(nc*nc*10) */
      /*    *\/ */
      /*   /\*      || squared_norm(localreaction) <=  (*norm_r* *norm_r/(nc*nc*1000))) *\/ */
      /*   /\*     && iter >=10) *\/ */
      /*   { */
      /*     /\* we  freeze the contact for n iterations*\/ */

      /*     freeze_contacts[contact] = options->iparam[SICONOS_NSGS_FREEZING_CONTACT]; */

      /*     DEBUG_EXPR( */
      /*         NV_display(localreaction, 3); NV_display(&reaction[contact * 3], 3); */
      /*         printf("light_error_2 = %e\n", light_error_2[contact]); */
      /*         printf("tmp_criteria1 = %e\n", tmp_criteria1); */
      /*         printf("tmp_criteria2 = %e\n", tmp_criteria2); */
      /*         printf("first criteria relative_convergence_criteria : light_error_2 <= " */
      /*                "tmp_criteria1 * squared_norm_localreaction ==> %e <= %e, bool =%i\n", */
      /*                light_error_2[contact], tmp_criteria1 * squared_norm_localreaction, */
      /*                relative_convergence_criteria); */
      /*         printf("second criteria :  squared_norm_localreaction <= tmp_criteria2 ==> %e " */
      /*                "<= %e, bool =%i \n", */
      /*                squared_norm_localreaction, tmp_criteria2, small_reaction_criteria); */
      /*         printf("Contact % i is freezed for %i steps\n", contact, */
      /*                options->iparam[SICONOS_NSGS_FREEZING_CONTACT]);); */
      /*   } */
      /* } */

      if (options->iparam[SICONOS_NSGS_FILTER_LOCAL_SOLUTION] ==
          SICONOS_NSGS_FILTER_LOCAL_SOLUTION_TRUE)
        acceptLocalReactionFiltered(localproblem_contact->dimension, local_opts_contact, contact, iter, reaction,
                                    localreaction);
      else
        acceptLocalReactionUnconditionally(contact, reaction, localreaction);
    }

    if (options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
        SICONOS_NSGS_ERROR_EVALUATION_LIGHT) {
      incr_error = calculateLightError(light_error_sum, nc, reaction, norm_r);
      hasNotConverged = nsgs_determine_convergence(incr_error, tolerance, iter, options);
    } else if (options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
               SICONOS_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
      incr_error = calculateLightError(light_error_sum, nc, reaction, norm_r);
      hasNotConverged = determine_convergence_with_full_final(
          problem, options, computeError, reaction, velocity, &tolerance, norm_q, incr_error,
          &full_error, iter);

      if (!(tolerance > 0.0)) {
        numerics_warning("cohesive_friction_3d_nsgs", "tolerance has to be positive!!");
        numerics_warning("cohesive_friction_3d_nsgs", "we stop the iterations");
        break;
      }

    } else if (options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
               SICONOS_NSGS_ERROR_EVALUATION_FULL) {
      full_error = calculateFullErrorAdaptiveInterval(problem, computeError, options, iter,
                                                      reaction, velocity, tolerance, norm_q);
      incr_error = full_error;
      hasNotConverged = nsgs_determine_convergence(full_error, tolerance, iter, options);
    }

    statsIterationCallback(problem, options, reaction, velocity, incr_error);
    if (verbose > 0) {
      frozen_contact = 0;
      if (options->iparam[SICONOS_NSGS_FREEZING_CONTACT] > 0) {
        for (unsigned int i = 0; i < nc; ++i) {
          if (freeze_contacts[i] > 0) {
            frozen_contact++;
          }
        }
      }
    }
    nsgs_print_iteration_stats(iter, incr_error, full_error, tolerance, SOLVER_TOL(options),
                               frozen_contact,  hasNotConverged, verbose);
  }
  free(light_error_2);

  /* Full criterium */
  if (options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
      SICONOS_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
    full_error = calculateFullErrorFinal(problem, options, computeError, reaction, velocity,
                                         tolerance, norm_q);

    hasNotConverged = nsgs_determine_convergence(full_error, SOLVER_TOL(options), iter, options);
  }
  /* if (iter == itermax) */
  /*   {             */
  /*     char fname[256];           */
  /*     fccounter++; */
  /*     snprintf(fname, sizeof(fname), "problem_%i_%i.dat", nc,itermax); */
  /*     FILE * file = fopen(fname, "w"); */
  /*     frictionContact_printInFile(problem, file); */
  /*     fclose(file); */
  /*   }  */

  /** return parameter values */
  /* SOLVER_TOL(options) = tolerance; */
  SET_SOLVER_RESIDUAL(options, full_error);
  SET_SOLVER_ITER_DONE(options, iter);

  /* Restore original local solver tolerance */
  tolerance_manager_restore_local(&tol_manager, local_opts_contact);

  /** Free memory **/
  if (problem->M->storageType == NM_SPARSE) {
    SBM_clear_block(problem->M->matrix1);
    SBM_clear(problem->M->matrix1);
    problem->M->matrix1 = matrix1;
  }
  localProblemFunctionToolkit->free_local_solver_contact(localproblem_contact, local_opts_contact);
  localProblemFunctionToolkit->free_local_solver_cohesion(localproblem_cohesion, local_opts_contact);

  cohesive_friction_3d_local_problem_free(localproblem_cohesion,  problem);


  localproblem_contact->M =  NULL;
  localproblem_contact->q =  NULL;
  localproblem_contact->mu =  NULL;
  frictionContactProblem_free(localproblem_contact);
  
  DEBUG_PRINTF("Final iteration: %d, error: %e\n", iter, error);
  DEBUG_END("cohesive_friction_3d_nsgs(...)\n");

  if (iter == itermax && hasNotConverged > 0) return NUMERICS_ERR_MAX_ITER;
  else
    {
      return (full_error <= tolerance) ? NUMERICS_OK : NUMERICS_ERR_DIVERGENCE;
    }    
}

/* ===========================================================================
 * Solver Registration
 * =========================================================================== */

static int cohesive_friction_3d_nsgs_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  if (!options) return NUMERICS_ERR_NULL_POINTER;

  return NUMERICS_OK;
}

static int cohesive_friction_3d_nsgs_solve_wrap(void* problem, double* reaction,
                                                double* velocity, SolverOptions* options) {
  return cohesive_friction_3d_nsgs((CohesiveFrictionContactProblem*)problem, reaction,
                                   velocity, options);
}

static void cohesive_friction_3d_nsgs_free_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  // Nothing to free specifically for NSGS
}

static void cohesive_friction_3d_nsgs_set_default(SolverOptions* options) {
  if (!options) return;

  options->iparam[SICONOS_NSGS_ERROR_EVALUATION_TYPE] =
      SICONOS_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL;
  options->iparam[SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE] =
      SICONOS_NSGS_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE;

  options->iparam[SICONOS_NSGS_SHUFFLE] = SICONOS_NSGS_SHUFFLE_FALSE;
  options->iparam[SICONOS_NSGS_SHUFFLE_SEED] = 0;
  options->iparam[SICONOS_NSGS_FREEZING_CONTACT] = 0;
  options->iparam[SICONOS_NSGS_FILTER_LOCAL_SOLUTION] =
      SICONOS_NSGS_FILTER_LOCAL_SOLUTION_FALSE;
  options->iparam[SICONOS_NSGS_RELAXATION] = SICONOS_NSGS_RELAXATION_FALSE;
  options->iparam[SICONOS_NSGS_ERROR_EVALUATION_FREQUENCY] = 0;

  SOLVER_TOL(options) = 1e-4;

  options->dparam[SICONOS_NSGS_INTERNAL_ERROR_RATIO] = 10.0;
  // Internal solver - allocate if needed
  if (options->numberOfInternalSolvers == 0) {
    options->numberOfInternalSolvers = 2;
    options->internalSolvers = calloc(2, sizeof(SolverOptions*));
  }
  assert(options->numberOfInternalSolvers == 2);
  options->internalSolvers[0] = solver_options_create(OC_NSN_GP_HYBRID);

  // to be changed for the cohesion solver.
  options->internalSolvers[1] = solver_options_create(SICONOS_COHESIVE_FRICTION_3D_PROJECTION);
}

/* Register the solver in the solver registry */
REGISTER_SOLVER(SICONOS_COHESIVE_FRICTION_3D_NSGS, "COHESIVE_FRICTION_3D_NSGS",
                "Non-smooth Gauss-Seidel for 3D Cohesive Friction Contact",
                cohesive_friction_3d_nsgs_init_wrap, cohesive_friction_3d_nsgs_solve_wrap,
                cohesive_friction_3d_nsgs_free_wrap,
                NULL, /* error function - use compute_error_velocity */
                cohesive_friction_3d_nsgs_set_default, 1000, /* default max iter */
                1e-4,                                        /* default tolerance */
                false /* not a local solver */);
