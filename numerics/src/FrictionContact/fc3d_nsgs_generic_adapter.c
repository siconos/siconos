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

/*!\file fc3d_nsgs_generic_adapter.c
 * \brief Adapter for fc3d NSGS using the generic NSGS framework
 *
 * This adapter allows fc3d problems to use the generic NSGS solver
 * by providing appropriate callbacks and initialization.
 */

#include "FrictionContactProblem.h"
#include "NumericsFwd.h"
#include "NumericsMatrix.h"
#include "fc3d_Solvers.h"
#include "nsgs_generic.h"
#include "Friction_cst.h"
#include "fc3d_local_problem_tools.h"
#include "fc3d_projection.h"
#include "fc3d_onecontact_nonsmooth_Newton_solvers.h"
#include "fc3d_Path.h"
#include "fc3d_NCPGlockerFixedPoint.h"
#include "fc3d_unitary_enumerative.h"
#include "numerics_verbose.h"
#include "fc3d_2NCP_Glocker.h"
#include "fc3d_compute_error.h"

/* Profiling support - define NSGS_PROFILE to enable */
//#define NSGS_PROFILE
#ifdef NSGS_PROFILE
#include <sys/time.h>
static double profile_update_time = 0.0;
static double profile_solve_time = 0.0;
static double profile_accept_time = 0.0;
static int profile_update_calls = 0;
static int profile_solve_calls = 0;
static int profile_accept_calls = 0;

static double profile_get_time(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

#define PROFILE_START() double __t_start = profile_get_time()
#define PROFILE_END(timer, counter) do { \
  timer += profile_get_time() - __t_start; \
  counter++; \
} while(0)

#define PROFILE_PRINT() do { \
  printf("\n=== NSGS Generic Profile ===\n"); \
  printf("update_local: %12.6f ms (%d calls, %.6f ms/call)\n", \
         profile_update_time * 1000, profile_update_calls, \
         profile_update_time * 1000 / (profile_update_calls > 0 ? profile_update_calls : 1)); \
  printf("solve_local:  %12.6f ms (%d calls, %.6f ms/call)\n", \
         profile_solve_time * 1000, profile_solve_calls, \
         profile_solve_time * 1000 / (profile_solve_calls > 0 ? profile_solve_calls : 1)); \
  printf("accept_local: %12.6f ms (%d calls, %.6f ms/call)\n", \
         profile_accept_time * 1000, profile_accept_calls, \
         profile_accept_time * 1000 / (profile_accept_calls > 0 ? profile_accept_calls : 1)); \
  printf("===========================\n\n"); \
} while(0)

#define PROFILE_RESET() do { \
  profile_update_time = profile_solve_time = profile_accept_time = 0.0; \
  profile_update_calls = profile_solve_calls = profile_accept_calls = 0; \
} while(0)

#else
#define PROFILE_START()
#define PROFILE_END(timer, counter)
#define PROFILE_PRINT()
#define PROFILE_RESET()
#endif

extern void fc3d_set_internalsolver_tolerance(FrictionContactProblem* problem,
                                              SolverOptions* options,
                                              SolverOptions* localsolver_options,
                                              double error);

/* Dimension-specific relaxation for 3D: z_new = omega*z_new + (1-omega)*z_old */
static void fc3d_relaxation_3(double* z_new, double* z_old, double omega) {
  z_new[0] = omega * z_new[0] + (1.0 - omega) * z_old[0];
  z_new[1] = omega * z_new[1] + (1.0 - omega) * z_old[1];
  z_new[2] = omega * z_new[2] + (1.0 - omega) * z_old[2];
}

/* Dimension-specific squared norm for 3D: ||z||^2 */
static double fc3d_squared_norm_3(double* z) {
  return z[0]*z[0] + z[1]*z[1] + z[2]*z[2];
}

/* Dimension-specific incremental error for 3D: ||z_new - z_old||^2 */
static double fc3d_incr_error_3(double* z_new, double* z_old) {
  double diff[3] = {z_new[0] - z_old[0], z_new[1] - z_old[1], z_new[2] - z_old[2]};
  return diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
}

/** Allocate and initialize local problem, return local solver function
 * 
 * \param[in] problem Global friction contact problem
 * \param[in] options NSGS solver options
 * \param[out] local_solver Output pointer for local solver function
 * \return Allocated local problem, or NULL on error
 */
static FrictionContactProblem* fc3d_nsgs_local_problem_new(FrictionContactProblem* problem, 
                                                            SolverOptions* options,
                                                            SolverPtr* local_solver) {
  FrictionContactProblem* localproblem = fc3d_local_problem_allocate(problem);
  if (!localproblem) {
    numerics_error("fc3d_nsgs_local_problem_new", "Failed to allocate local problem");
    return NULL;
  }
  
  SolverOptions* localsolver_options = options->internalSolvers[0];
  SolverPtr solver = NULL;
  
  switch (localsolver_options->solverId) {
    case SICONOS_ONECONE_ProjectionOnConeWithDiagonalization:
      solver = &fc3d_projectionWithDiagonalization_solve;
      fc3d_projection_initialize(problem, localproblem);
      break;
    case SICONOS_ONECONE_ProjectionOnCone:
      solver = &fc3d_projectionOnCone_solve;
      fc3d_projection_initialize(problem, localproblem);
      break;
    case SICONOS_ONECONE_ProjectionOnConeWithLocalIteration:
      solver = &fc3d_projectionOnConeWithLocalIteration_solve;
      fc3d_projectionOnConeWithLocalIteration_initialize(problem, localproblem, localsolver_options);
      break;
    case SICONOS_ONECONE_ProjectionOnConeWithRegularization:
      solver = &fc3d_projectionOnCone_solve;
      fc3d_projection_initialize_with_regularization(problem, localproblem);
      break;
    case SICONOS_ONECONE_NSN:
      solver = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
      break;
    case SICONOS_ONECONE_NSN_GP:
      solver = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
      break;
    case SICONOS_ONECONE_NSN_GP_HYBRID:
      solver = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
      break;
    case SICONOS_FRICTION_3D_NCPGlockerFBNewton:
      solver = &fc3d_onecontact_nonsmooth_Newton_solvers_solve;
      fc3d_onecontact_nonsmooth_Newton_solvers_initialize(problem, localproblem, localsolver_options);
      break;
    case SICONOS_FRICTION_3D_NCPGlockerFBPATH:
      solver = &fc3d_Path_solve;
      fc3d_Path_initialize(problem, localproblem, localsolver_options);
      break;
    case SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint:
      solver = &fc3d_FixedP_solve;
      fc3d_FixedP_initialize(problem, localproblem, localsolver_options);
      break;
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder:
      solver = &fc3d_projectionOnCylinder_solve;
      fc3d_projectionOnCylinder_initialize(problem, localproblem, options);
      break;
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration:
      solver = &fc3d_projectionOnCylinderWithLocalIteration_solve;
      fc3d_projectionOnCylinderWithLocalIteration_initialize(problem, localproblem, options, localsolver_options);
      break;
    case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC:
      solver = &fc3d_unitary_enumerative_solve;
      fc3d_unitary_enumerative_initialize(localproblem);
      break;
    case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU:
      solver = &fc3d_unitary_enumerative_solve;
      fc3d_unitary_enumerative_initialize(localproblem);
      break;
    default:
      fc3d_local_problem_free(localproblem, problem);
      numerics_error("fc3d_nsgs_local_problem_new", "Unknown internal solver: %s",
                     solver_options_id_to_name(localsolver_options->solverId));
      return NULL;
  }
  
  *local_solver = solver;
  return localproblem;
}

/** Free local problem */
static void fc3d_nsgs_local_problem_free(FrictionContactProblem* localproblem, 
                                          FrictionContactProblem* problem) {
  if (localproblem) {
    fc3d_local_problem_free(localproblem, problem);
  }
}

/** Update local problem wrapper - maps generic NSGS to fc3d */
static void fc3d_nsgs_update_wrapper(unsigned int block, void* problem, void* local_problem,
                                     double* var_z, SolverOptions* options) {
  if (!local_problem) {
    numerics_error("fc3d_nsgs_update_wrapper", "Local problem is NULL");
    return;
  }
  PROFILE_START();
  fc3d_nsgs_update(block, (FrictionContactProblem*)problem,
                   (FrictionContactProblem*)local_problem, var_z, options);
  PROFILE_END(profile_update_time, profile_update_calls);
}

/** Solve local problem wrapper */
static void fc3d_nsgs_solve_local_wrapper(void* local_problem, SolverOptions* localsolver_options,
                                          double* var_z_local, double* localsolver_options_data) {
  (void)localsolver_options_data;
  
  SolverPtr local_solver = (SolverPtr)localsolver_options->solverData;
  if (local_solver && local_problem) {
    PROFILE_START();
    local_solver((FrictionContactProblem*)local_problem, var_z_local, localsolver_options);
    PROFILE_END(profile_solve_time, profile_solve_calls);
  } else {
    numerics_error("fc3d_nsgs_solve_local_wrapper", "Local solver or local problem not initialized");
  }
}

/** Compute error wrapper - maps generic to fc3d */
static double fc3d_nsgs_compute_error_wrapper(void* problem, double* var_z,
                                              double* var_x, SolverOptions* options) {
  double error = 0.0;
  double norm_q = 0.0;
  fc3d_compute_error((FrictionContactProblem*)problem, var_z, var_x, 
                     options->dparam[SICONOS_DPARAM_TOL], options, norm_q, &error);
  return error;
}

/** Set local tolerance wrapper */
static void fc3d_nsgs_set_tol_wrapper(void* problem, SolverOptions* options,
                                      SolverOptions* localsolver_options, double error) {
  fc3d_set_internalsolver_tolerance((FrictionContactProblem*)problem, options,
                                    localsolver_options, error);
}

/** Accept local solution wrapper with filtering */
static void fc3d_nsgs_accept_local_wrapper(void* local_problem, SolverOptions* options,
                                           unsigned int block, int iter,
                                           double* var_z_global, double* var_z_local) {
  (void)local_problem;

  PROFILE_START();
  double local_residual = options->dparam[SICONOS_DPARAM_RESIDU];
  
  /* DEBUG: Always print for first few iterations */
  if (iter <= 2 && block < 3) {
    printf("DEBUG accept: block=%d iter=%d local=[%.4e, %.4e, %.4e] residual=%.4e ",
           block, iter, var_z_local[0], var_z_local[1], var_z_local[2], local_residual);
  }
  
  if (isnan(local_residual) || isinf(local_residual) || local_residual > 1.0) {
    numerics_printf("Discard local reaction for block %i at iteration %i with local_error = %e",
                    block, iter, local_residual);
    if (iter <= 2 && block < 3) printf("[DISCARDED]\n");
    PROFILE_END(profile_accept_time, profile_accept_calls);
    return;
  }

  var_z_global[block * 3 + 0] = var_z_local[0];
  var_z_global[block * 3 + 1] = var_z_local[1];
  var_z_global[block * 3 + 2] = var_z_local[2];
  
  if (iter <= 2 && block < 3) {
    printf("[ACCEPTED] global[%d]=[%.4e, %.4e, %.4e]\n",
           block, var_z_global[block*3], var_z_global[block*3+1], var_z_global[block*3+2]);
  }
  
  PROFILE_END(profile_accept_time, profile_accept_calls);
}

/** Main solver using generic NSGS framework */
void fc3d_nsgs_generic(FrictionContactProblem* problem, double* reaction, double* velocity,
                       int* info, SolverOptions* options) {
  if (!problem || !reaction || !info || !options) {
    numerics_error("fc3d_nsgs_generic", "Invalid input arguments");
    return;
  }
  
  if (options->numberOfInternalSolvers < 1) {
    numerics_error("fc3d_nsgs_generic",
                   "The NSGS method needs options for the internal solvers");
    *info = 1;
    return;
  }

  /* Extract diagonal blocks for efficient access (for NM_SPARSE matrices) */
  void* original_matrix1 = nsgs_generic_extract_diagonal_blocks(problem->M, problem->dimension);

  /* Create local problem and get local solver function */
  SolverPtr local_solver = NULL;
  FrictionContactProblem* localproblem = fc3d_nsgs_local_problem_new(problem, options, &local_solver);
  if (!localproblem) {
    *info = 1;
    return;
  }
  
  /* Store local_solver in internal solver options for the wrapper to access */
  SolverOptions* localsolver_options = options->internalSolvers[0];
  localsolver_options->solverData = (void*)local_solver;

  /* Setup problem data */
  NSGSProblemData problem_data = {
    .nb_blocks = problem->numberOfContacts,
    .q = problem->q,
    .M = problem->M,
    .mu = problem->mu,
    .mu_r = NULL,
    .storage_type = problem->M->storageType,
    .dimension = 3
  };

  /* Dimension-specific copy function for 3D */
  NSGSCopyLocal copy_local_3 = NULL;  /* Use default BLAS copy in nsgs_solve */
  
  /* Setup generic NSGS toolkit */
  NSGSLocalToolkit toolkit = {
    .update_local_problem = fc3d_nsgs_update_wrapper,
    .copy_local = copy_local_3,
    .solve_local = fc3d_nsgs_solve_local_wrapper,
    .compute_error = fc3d_nsgs_compute_error_wrapper,
    .incremental_error = fc3d_incr_error_3,
    .accept_local = fc3d_nsgs_accept_local_wrapper,
    .check_convergence = NULL,
    .alloc_local = NULL,
    .set_local_tol = fc3d_nsgs_set_tol_wrapper,
    .stats_callback = NULL,
    .relaxation = fc3d_relaxation_3,
    .squared_norm = fc3d_squared_norm_3,
    .should_freeze = NULL,
    .alloc_freezing = NULL,  /* Use default allocation in nsgs_solve */
    .alloc_shuffled = NULL,  /* Use default allocation in nsgs_solve */
    .localproblem = localproblem,
    .verbose = 2,                    /* Print every iteration */
    .user_tolerance = options->dparam[SICONOS_DPARAM_TOL],  /* Store user tolerance */
    .error_eval_freq = 0,            /* Full error every iteration (0=always) */
    .dimension = 3,
    .omega = options->dparam[SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE],
    .use_freezing = options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0,
    .use_shuffling = options->iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] != SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE,
    .use_relaxation = options->iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION] == SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE,
    .use_incremental_error = 1,
    .filter_local_sol = options->iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION],
    .error_eval_type = options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION],
    .freezing_iter = options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT]
  };

  /* Reset profiling counters */
  PROFILE_RESET();
  
  /* Call generic solver */
  nsgs_solve(problem, reaction, velocity, info, options, &toolkit, &problem_data);

  /* Print profiling results */
  PROFILE_PRINT();

  /* Cleanup */
  fc3d_nsgs_local_problem_free(localproblem, problem);
  localsolver_options->solverData = NULL;

  /* Free extracted diagonal blocks and restore original matrix1 */
  nsgs_generic_free_diagonal_blocks(problem->M, original_matrix1);
}
