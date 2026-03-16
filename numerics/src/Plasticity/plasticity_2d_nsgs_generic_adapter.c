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

/*!\file plasticity_2d_nsgs_generic_adapter.c
 * \brief Generic NSGS adapter for 2D plasticity problems
 *
 * This adapter allows plasticity_2d problems to use the generic NSGS solver
 * by providing appropriate callbacks and initialization.
 */

#include "PlasticityProblem.h"
#include "NumericsFwd.h"
#include "NumericsMatrix.h"
#include "plasticity_2d_solvers.h"
#include "nsgs_generic.h"
#include "Plasticity_options.h"
#include "plasticity_2d_local_problem_tools.h"
#include "plasticity_2d_projection.h"
#include "plasticity_2d_onecone_nonsmooth_Newton_solvers.h"
#include "plasticity_2d_compute_error.h"
#include "numerics_verbose.h"
#include "numerics_errors.h"

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

extern int plasticity_2d_set_internalsolver_tolerance(PlasticityProblem* problem,
                                                       SolverOptions* options,
                                                       SolverOptions* localsolver_options,
                                                       double error);

/* Forward declaration for the update function */
static int  plasticity_2d_nsgs_update(int cone, PlasticityProblem* problem,
                                      PlasticityProblem* localproblem, double* stress,
                                      SolverOptions* options);

/* Dimension-specific relaxation for 3D: z_new = omega*z_new + (1-omega)*z_old */
static void plasticity_2d_relaxation_3(double* z_new, double* z_old, double omega) {
  z_new[0] = omega * z_new[0] + (1.0 - omega) * z_old[0];
  z_new[1] = omega * z_new[1] + (1.0 - omega) * z_old[1];
  z_new[2] = omega * z_new[2] + (1.0 - omega) * z_old[2];
}

/* Dimension-specific squared norm for 3D: ||z||^2 */
static double plasticity_2d_squared_norm_3(double* z) {
  return z[0]*z[0] + z[1]*z[1] + z[2]*z[2];
}

/* Dimension-specific incremental error for 3D: ||z_new - z_old||^2 */
static double plasticity_2d_incr_error_3(double* z_new, double* z_old) {
  double diff[3] = {z_new[0] - z_old[0], z_new[1] - z_old[1], z_new[2] - z_old[2]};
  return diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
}

/** Allocate and initialize local problem, return local solver function
 * 
 * \param[in] problem Global plasticity problem
 * \param[in] options NSGS solver options
 * \param[out] local_solver Output pointer for local solver function
 * \return Allocated local problem, or NULL on error
 */
static PlasticityProblem* plasticity_2d_nsgs_local_problem_new(PlasticityProblem* problem, 
                                                                SolverOptions* options,
                                                                SolverPtr* local_solver) {
  PlasticityProblem* localproblem = plasticity_2d_local_problem_allocate(problem);
  if (!localproblem) return NULL;
  
  SolverOptions* localsolver_options = options->internalSolvers[0];
  SolverPtr solver = NULL;
  
  /* Initialize based on model type */
  if (problem->model_type == PLASTICITY_MODEL_DRUCKER_PRAGER) {
    switch (localsolver_options->solverId) {
      case PLASTICITY_2D_ONECONE_ProjectionOnCone:
        solver = &plasticity_2d_projectionOnCone_solve;
        plasticity_2d_projection_initialize(problem, localproblem);
        /* Allocate dWork for rho values (needed by projection solver) */
        if (!localsolver_options->dWork || localsolver_options->dWorkSize < problem->numberOfCones) {
          localsolver_options->dWork = (double*)realloc(localsolver_options->dWork, 
                                                         problem->numberOfCones * sizeof(double));
          localsolver_options->dWorkSize = problem->numberOfCones;
        }
        for (size_t i = 0; i < (size_t)problem->numberOfCones; i++) {
          localsolver_options->dWork[i] = 1.0;
        }
        break;
      case PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration:
        solver = &plasticity_2d_projectionOnConeWithLocalIteration_solve;
        plasticity_2d_projectionOnConeWithLocalIteration_initialize(problem, localproblem,
                                                                    localsolver_options);
        break;
      case PLASTICITY_2D_ONECONE_NSN:
      case PLASTICITY_2D_ONECONE_NSN_GP:
      case PLASTICITY_2D_ONECONE_NSN_GP_HYBRID:
        solver = &plasticity_2d_onecone_nonsmooth_Newton_solvers_solve;
        plasticity_2d_onecone_nonsmooth_Newton_solvers_initialize(problem, localproblem,
                                                                  localsolver_options);
        break;
      default:
        plasticity_2d_local_problem_free(localproblem, problem);
        int error  = numerics_error("plasticity_2d_nsgs_local_problem_new", 
                       "Unknown internal solver: %s",
                       solver_options_id_to_name(localsolver_options->solverId));
        return NULL;
    }
  } else if (problem->model_type == PLASTICITY_MODEL_VON_MISES) {
    switch (localsolver_options->solverId) {
      case PLASTICITY_2D_ONECONE_ProjectionOnCone:
        solver = &plasticity_2d_projectionOnVonMises_solve;
        plasticity_2d_projection_initialize(problem, localproblem);
        break;
      case PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration:
        solver = &plasticity_2d_projectionOnVonMises_solve;
        plasticity_2d_projectionOnConeWithLocalIteration_initialize(problem, localproblem,
                                                                    localsolver_options);
        break;
      default:
        plasticity_2d_local_problem_free(localproblem, problem);
        int error = numerics_error("plasticity_2d_nsgs_local_problem_new", 
                       "Unknown internal solver for Von Mises: %s",
                       solver_options_id_to_name(localsolver_options->solverId));
        return NULL;
    }
  } else {
    plasticity_2d_local_problem_free(localproblem, problem);
    int error = numerics_error("plasticity_2d_nsgs_local_problem_new", 
                   "Unknown plasticity model type");
    return NULL;
  }
  
  *local_solver = solver;
  return localproblem;
}

/** Free local problem */
static void plasticity_2d_nsgs_local_problem_free(PlasticityProblem* localproblem, 
                                                   PlasticityProblem* problem) {
  if (localproblem) {
    plasticity_2d_local_problem_free(localproblem, problem);
  }
}

/** Update local problem wrapper - maps generic NSGS to plasticity_2d */
static int plasticity_2d_nsgs_update_wrapper(unsigned int block, void* problem, void* local_problem,
                                               double* var_z, SolverOptions* options) {
  if (!local_problem) {
    return numerics_error("plasticity_2d_nsgs_update_wrapper", "Local problem is NULL");
  }
  PROFILE_START();
  int info = plasticity_2d_nsgs_update(block, (PlasticityProblem*)problem,
                            (PlasticityProblem*)local_problem, var_z, options);
  PROFILE_END(profile_update_time, profile_update_calls);
  return info;
}

/** Solve local problem wrapper */
static int  plasticity_2d_nsgs_solve_local_wrapper(void* local_problem, 
                                                    SolverOptions* localsolver_options,
                                                    double* var_z_local, 
                                                    double* localsolver_options_data) {
  (void)localsolver_options_data;
  
  SolverPtr local_solver = (SolverPtr)localsolver_options->solverData;
  if (local_solver && local_problem) {
    PROFILE_START();
    local_solver((PlasticityProblem*)local_problem, var_z_local, localsolver_options);
    /* For solvers that don't set residual, set it to a valid value (0) */
    if (isnan(localsolver_options->dparam[SICONOS_DPARAM_RESIDU]) ||
        isinf(localsolver_options->dparam[SICONOS_DPARAM_RESIDU])) {
      localsolver_options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;
    }
    PROFILE_END(profile_solve_time, profile_solve_calls);
  } else {
    return numerics_error("plasticity_2d_nsgs_solve_local_wrapper", 
                   "Local solver or local problem not initialized");
  }
  return 0;
}

/** Compute error wrapper - maps generic to plasticity_2d */
static double plasticity_2d_nsgs_compute_error_wrapper(void* problem, double* var_z,
                                                        double* var_x, SolverOptions* options) {
  double error = 0.0;
  double norm_q = 0.0;
  plasticity_2d_compute_error((PlasticityProblem*)problem, var_z, var_x, 
                               options->dparam[SICONOS_DPARAM_TOL], options, norm_q, &error);
  return error;
}

/** Set local tolerance wrapper */
static void plasticity_2d_nsgs_set_tol_wrapper(void* problem, SolverOptions* options,
                                                SolverOptions* localsolver_options, 
                                                double error) {
  plasticity_2d_set_internalsolver_tolerance((PlasticityProblem*)problem, options,
                                             localsolver_options, error);
}

/** Accept local solution wrapper with filtering */
static int  plasticity_2d_nsgs_accept_local_wrapper(void* local_problem, 
                                                     SolverOptions* options,
                                                     unsigned int block, int iter,
                                                     double* var_z_global, 
                                                     double* var_z_local) {
  (void)local_problem;
  (void)iter;

  PROFILE_START();
  double local_residual = options->dparam[SICONOS_DPARAM_RESIDU];
  if (isnan(local_residual) || isinf(local_residual) || local_residual > 1.0) {
    numerics_printf("Discard local stress for block %i at iteration %i with local_error = %e",
                    block, iter, local_residual);
    PROFILE_END(profile_accept_time, profile_accept_calls);
    return 0;
  }

  var_z_global[block * 3 + 0] = var_z_local[0];
  var_z_global[block * 3 + 1] = var_z_local[1];
  var_z_global[block * 3 + 2] = var_z_local[2];
  PROFILE_END(profile_accept_time, profile_accept_calls);
  return 0;  
}

/** Update local problem for a specific cone/block
 * 
 * This is the same implementation as in plasticity_2d_nsgs.c
 */
static int plasticity_2d_nsgs_update(int cone, PlasticityProblem* problem,
				     PlasticityProblem* localproblem, double* stress,
				     SolverOptions* options) {
  (void)options;
  
  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  plasticity_2d_local_problem_fill_M(problem, localproblem, cone);

  /* Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products
     MLocal.stressBlock, excluding the block corresponding to the current cone. */
  plasticity_2d_local_problem_compute_q(problem, localproblem, stress, cone);

  /* coefficient for current block - handle both model types */
  if (problem->model_type == PLASTICITY_MODEL_DRUCKER_PRAGER) {
    localproblem->model.drucker_prager->eta[0] = problem->model.drucker_prager->eta[cone];
    localproblem->model.drucker_prager->theta[0] = problem->model.drucker_prager->theta[cone];
  } else if (problem->model_type == PLASTICITY_MODEL_VON_MISES) {
    localproblem->model.von_mises->sigma_y[0] = problem->model.von_mises->sigma_y[cone];
  }
  return 0;  
}

/** Main solver using generic NSGS framework */
void plasticity_2d_nsgs_generic(PlasticityProblem* problem, double* stress, 
                                 double* plastic_strain_rate, int* info, 
                                 SolverOptions* options) {
  if (!problem || !stress || !info || !options) {
    *info = numerics_error("plasticity_2d_nsgs_generic", "Invalid input arguments");
    return;
  }
  
  if (options->numberOfInternalSolvers < 1) {
    *info = numerics_error("plasticity_2d_nsgs_generic",
                   "The NSGS method needs options for the internal solvers");
    return;
  }

  /* Extract diagonal blocks for efficient access (for NM_SPARSE matrices) */
  void* original_matrix1 = nsgs_generic_extract_diagonal_blocks(problem->M, problem->dimension);

  /* Create local problem and get local solver function */
  SolverPtr local_solver = NULL;
  PlasticityProblem* localproblem = plasticity_2d_nsgs_local_problem_new(problem, options, 
                                                                          &local_solver);
  if (!localproblem) {
    *info = 1;
    return;
  }
  
  /* Store local_solver in internal solver options for the wrapper to access */
  SolverOptions* localsolver_options = options->internalSolvers[0];
  localsolver_options->solverData = (void*)local_solver;
  


  /* Setup problem data - note: plasticity doesn't use mu/mu_r like friction contact */
  NSGSProblemData problem_data = {
    .nb_blocks = problem->numberOfCones,
    .q = problem->q,
    .M = problem->M,
    .mu = NULL,  /* Not used for plasticity */
    .mu_r = NULL,
    .storage_type = problem->M->storageType,
    .dimension = 3
  };

  /* Setup generic NSGS toolkit */
  NSGSLocalToolkit toolkit = {
    .update_local_problem = plasticity_2d_nsgs_update_wrapper,
    .copy_local = NULL,  /* Use default BLAS copy in nsgs_solve */
    .solve_local = plasticity_2d_nsgs_solve_local_wrapper,
    .compute_error = plasticity_2d_nsgs_compute_error_wrapper,
    .incremental_error = plasticity_2d_incr_error_3,
    .accept_local = plasticity_2d_nsgs_accept_local_wrapper,
    .check_convergence = NULL,
    .alloc_local = NULL,
    .set_local_tol = plasticity_2d_nsgs_set_tol_wrapper,
    .stats_callback = NULL,
    .relaxation = plasticity_2d_relaxation_3,
    .squared_norm = plasticity_2d_squared_norm_3,
    .should_freeze = NULL,
    .alloc_freezing = NULL,  /* Use default allocation in nsgs_solve */
    .alloc_shuffled = NULL,  /* Use default allocation in nsgs_solve */
    .localproblem = localproblem,
    .verbose = 2,                    /* Print every iteration */
    .user_tolerance = options->dparam[SICONOS_DPARAM_TOL],  /* Store user tolerance */
    .error_eval_freq = 0,            /* Full error every iteration (0=always) */
    .dimension = 3,
    .omega = options->dparam[PLASTICITY_NSGS_RELAXATION_VALUE],
    .use_freezing = options->iparam[PLASTICITY_NSGS_FREEZING_CONE] > 0,
    .use_shuffling = options->iparam[PLASTICITY_NSGS_SHUFFLE] != PLASTICITY_NSGS_SHUFFLE_FALSE,
    .use_incremental_error = 1,
    .use_relaxation = options->iparam[PLASTICITY_NSGS_RELAXATION] == PLASTICITY_NSGS_RELAXATION_TRUE,
    .filter_local_sol = options->iparam[PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION],
    .error_eval_type = options->iparam[PLASTICITY_IPARAM_ERROR_EVALUATION],
    .freezing_iter = options->iparam[PLASTICITY_NSGS_FREEZING_CONE]
  };

  /* Reset profiling counters */
  PROFILE_RESET();
  
  /* Call generic solver */
  nsgs_solve(problem, stress, plastic_strain_rate, info, options, &toolkit, &problem_data);

  /* Print profiling results */
  PROFILE_PRINT();

  /* Cleanup */
  plasticity_2d_nsgs_local_problem_free(localproblem, problem);
  localsolver_options->solverData = NULL;

  /* Free extracted diagonal blocks and restore original matrix1 */
  nsgs_generic_free_diagonal_blocks(problem->M, original_matrix1);
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers PLASTICITY_2D_NSGS_GENERIC in the global solver registry
 */

#include "solver_registry.h"

static void plasticity_2d_nsgs_generic_set_default(SolverOptions* options) {
  /* Delegate to the regular NSGS set_default since options are the same */
  extern void plasticity_2d_nsgs_set_default(SolverOptions*);
  plasticity_2d_nsgs_set_default(options);
}

static int plasticity_2d_nsgs_generic_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int plasticity_2d_nsgs_generic_solve_wrap(void* problem, double* reaction,
                                                  double* velocity, SolverOptions* options) {
  int info = NUMERICS_OK;
  plasticity_2d_nsgs_generic((PlasticityProblem*)problem, reaction, velocity, &info, options);
  return info;
}

static void plasticity_2d_nsgs_generic_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(PLASTICITY_2D_NSGS_GENERIC, "PLASTICITY_2D_NSGS_GENERIC",
                "Non-smooth Gauss-Seidel for 2D Plasticity (generic framework)",
                plasticity_2d_nsgs_generic_init_wrap,
                plasticity_2d_nsgs_generic_solve_wrap,
                plasticity_2d_nsgs_generic_free_wrap,
                NULL,  /* error function */
                plasticity_2d_nsgs_generic_set_default,  /* set_default */
                1000,  /* default_max_iter */
                1e-4,  /* default_tol */
                0      /* is_local_solver */);
