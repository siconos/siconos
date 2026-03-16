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

#include "PlasticityProblem.h"  // for PlasticityProblem
#include "NumericsArrays.h"        // for uint_shuffle
#include "NumericsFwd.h"           // for SolverOptions
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "Plasticity_options.h"  // for SICONOS_FRICTI...
#include "SiconosBlas.h"     // for cblas_dnrm2
#include "SolverOptions.h"   // for SolverOptions
#include "SparseBlockMatrix.h"
#include "plasticity_2d_compute_error.h"                     // for fc3d_compute_e..
#include "plasticity_2d_local_problem_tools.h"               // for fc3d_local_pro..
#include "plasticity_2d_onecone_nonsmooth_Newton_solvers.h"  //
#include "plasticity_2d_projection.h"                        // for fc3d_projectio...
#include "plasticity_2d_solvers.h"
#include "numerics_verbose.h"

/* Solver registration system */
#include "solver_registry.h"
#include "numerics_errors.h"  // for numerics_printf
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */

#include "op3x3.h"
#include "siconos_debug.h"  // for DEBUG_EXPR

// #define FCLIB_OUTPUT

#ifdef FCLIB_OUTPUT
static int fccounter = -1;
#include "fclib_interface.h"
#endif

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

/* static void fake_compute_error_nsgs(PlasticityProblem* problem, double *stress,
 * double *strainrate, double tolerance, SolverOptions  *options,  double* error) */
/* { */
/*   int n = 3 * problem->numberOfCones; */
/*   *error = 0.; */
/*   int i, m; */
/*   m = 5 * n / 3; */
/*   double err = INFINITY; */
/*   for (i = 0 ; i < m ; ++i) */
/*   { */
/*     *error += Compute_NCP_error1(i, err); */
/*   } */
/* } */

static inline void performRelaxation_3(double localstress[3], double *oldstress,
                                       double omega) {
  localstress[0] = omega * localstress[0] + (1.0 - omega) * oldstress[0];
  localstress[1] = omega * localstress[1] + (1.0 - omega) * oldstress[1];
  localstress[2] = omega * localstress[2] + (1.0 - omega) * oldstress[2];
}

static inline double light_error_squared_3(double localstress[3], double *oldstress) {
  return (pow(oldstress[0] - localstress[0], 2) +
          pow(oldstress[1] - localstress[1], 2) +
          pow(oldstress[2] - localstress[2], 2));
}

static inline double squared_norm_3(double localstress[3]) {
  return (pow(localstress[0], 2) + pow(localstress[1], 2) + pow(localstress[2], 2));
}

static double calculateLightError(double light_error_sum, unsigned int nc, double *stress,
                                  double *norm_r) {
  double error = sqrt(light_error_sum);
  *norm_r = cblas_dnrm2(nc * 3, stress, 1);
  if (fabs(*norm_r) > DBL_EPSILON) error /= (*norm_r);
  return error;
}

static void acceptLocalReactionUnconditionally(unsigned int cone, double *stress,
                                               double localstress[3]) {
  memcpy(&stress[cone * 3], localstress, sizeof(double) * 3);
}

static void statsIterationCallback(PlasticityProblem *problem, SolverOptions *options,
                                   double *stress, double *strainrate, double error) {
  if (options->callback) {
    options->callback->collectStatsIteration(
        options->callback->env, problem->numberOfCones * 3, stress, strainrate, error, NULL);
  }
}

static void plasticity_2d_nsgs_update(int cone, PlasticityProblem *problem,
                             PlasticityProblem *localproblem, double *stress,
                             SolverOptions *options) {
  /* Build a local problem for a specific cone
     stress corresponds to the global vector (size n) of the global problem.
  */
  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific cone
     stress corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  plasticity_2d_local_problem_fill_M(problem, localproblem, cone);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products
     MLocal.stressBlock, excluding the block corresponding to the current cone. ****/
  plasticity_2d_local_problem_compute_q(problem, localproblem, stress, cone);

  /* coefficient for current block - handle both model types */
  if (problem->model_type == PLASTICITY_MODEL_DRUCKER_PRAGER) {
    localproblem->model.drucker_prager->eta[0] = problem->model.drucker_prager->eta[cone];
    localproblem->model.drucker_prager->theta[0] = problem->model.drucker_prager->theta[cone];
  } else if (problem->model_type == PLASTICITY_MODEL_VON_MISES) {
    localproblem->model.von_mises->sigma_y[0] = problem->model.von_mises->sigma_y[cone];
  }
}

int plasticity_2d_nsgs_initialize_local_solver(
    struct LocalPLASTICITY_2DProblemFunctionToolkit *local_function_toolkit,
    plasticity_2d_ComputeErrorPtr *computeError, PlasticityProblem *problem,
    PlasticityProblem *localproblem, SolverOptions *options) {
  SolverOptions *localsolver_options = options->internalSolvers[0];

  *computeError = (plasticity_2d_ComputeErrorPtr)&plasticity_2d_compute_error;

  if (problem->dimension == 3) {
    local_function_toolkit->copy_local_reaction = cpy3;
    local_function_toolkit->perform_relaxation = &performRelaxation_3;
    local_function_toolkit->light_error_squared = &light_error_squared_3;
    local_function_toolkit->squared_norm = &squared_norm_3;
  }
  if (problem->model_type == PLASTICITY_MODEL_DRUCKER_PRAGER) {
    /** Connect to local solver */
    switch (localsolver_options->solverId) {
      case PLASTICITY_2D_ONECONE_ProjectionOnCone: {
        local_function_toolkit->local_solver = &plasticity_2d_projectionOnCone_solve;
        local_function_toolkit->update_local_problem = &plasticity_2d_nsgs_update;
        local_function_toolkit->free_local_solver = &plasticity_2d_projection_free;
        plasticity_2d_projection_initialize(problem, localproblem);
        break;
      }
      case PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration: {
        local_function_toolkit->local_solver =
            &plasticity_2d_projectionOnConeWithLocalIteration_solve;
        local_function_toolkit->update_local_problem = &plasticity_2d_nsgs_update;
        local_function_toolkit->free_local_solver =
            &plasticity_2d_projectionOnConeWithLocalIteration_free;
        plasticity_2d_projectionOnConeWithLocalIteration_initialize(problem, localproblem,
                                                                    localsolver_options);
        break;
      }
      /* Newton solver (Alart-Curnier) */
      case PLASTICITY_2D_ONECONE_NSN: {
        local_function_toolkit->local_solver =
            &plasticity_2d_onecone_nonsmooth_Newton_solvers_solve;
        local_function_toolkit->update_local_problem =
            &plasticity_2d_onecone_nonsmooth_Newton_AC_update;
        local_function_toolkit->free_local_solver =
            &plasticity_2d_onecone_nonsmooth_Newton_solvers_free;
        plasticity_2d_onecone_nonsmooth_Newton_solvers_initialize(problem, localproblem,
                                                                  localsolver_options);
        break;
      }
      case PLASTICITY_2D_ONECONE_NSN_GP: {
        local_function_toolkit->local_solver =
            &plasticity_2d_onecone_nonsmooth_Newton_solvers_solve;
        local_function_toolkit->update_local_problem =
            &plasticity_2d_onecone_nonsmooth_Newton_AC_update;
        local_function_toolkit->free_local_solver =
            &plasticity_2d_onecone_nonsmooth_Newton_solvers_free;
        plasticity_2d_onecone_nonsmooth_Newton_solvers_initialize(problem, localproblem,
                                                                  localsolver_options);
        break;
      }
      case PLASTICITY_2D_ONECONE_NSN_GP_HYBRID: {
        local_function_toolkit->local_solver =
            &plasticity_2d_onecone_nonsmooth_Newton_solvers_solve;
        local_function_toolkit->update_local_problem =
            &plasticity_2d_onecone_nonsmooth_Newton_AC_update;
        local_function_toolkit->free_local_solver =
            &plasticity_2d_onecone_nonsmooth_Newton_solvers_free;
        plasticity_2d_onecone_nonsmooth_Newton_solvers_initialize(problem, localproblem,
                                                                  localsolver_options);
        break;
      }
      default: {
        return numerics_error("plasticity_2d_nsgs_initialize_local_solver",
                       "Numerics, plasticity_2d_nsgs failed. Unknown internal solver : %s.\n",
                       solver_options_id_to_name(localsolver_options->solverId));
      }
    }
  } else if (problem->model_type == PLASTICITY_MODEL_VON_MISES) {
    switch (localsolver_options->solverId) {
      case PLASTICITY_2D_ONECONE_ProjectionOnCone: {
        local_function_toolkit->local_solver = &plasticity_2d_projectionOnVonMises_solve;
        local_function_toolkit->update_local_problem = &plasticity_2d_nsgs_update;
        local_function_toolkit->free_local_solver = &plasticity_2d_projection_free;
        plasticity_2d_projection_initialize(problem, localproblem);
        break;
      }
      case PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration: {
        local_function_toolkit->local_solver = &plasticity_2d_projectionOnVonMises_solve;
        local_function_toolkit->update_local_problem = &plasticity_2d_nsgs_update;
        local_function_toolkit->free_local_solver =
            &plasticity_2d_projectionOnConeWithLocalIteration_free;
        plasticity_2d_projectionOnConeWithLocalIteration_initialize(problem, localproblem,
                                                                    localsolver_options);
        break;
      }
      default: {
        numerics_printf_verbose(1,"plasticity_2d_nsgs_initialize_local_solver, Numerics, plasticity_2d_nsgs failed. Unknown internal solver : %s.\n",
                       solver_options_id_to_name(localsolver_options->solverId));
      }
    }
  } else {
    return numerics_error("plasticity_2d_nsgs_initialize_local_solver",
		   "Numerics, plasticity_2d_nsgs failed. Unknown plasticity model \n");
    }

    return 0;
    }

static unsigned int *allocShuffledCones(PlasticityProblem *problem,
                                        SolverOptions *options) {
  unsigned int *scones = 0;
  unsigned int nc = problem->numberOfCones;
  if (options->iparam[PLASTICITY_NSGS_SHUFFLE] == PLASTICITY_NSGS_SHUFFLE_TRUE ||
      options->iparam[PLASTICITY_NSGS_SHUFFLE] == PLASTICITY_NSGS_SHUFFLE_TRUE_EACH_LOOP) {
    if (options->iparam[PLASTICITY_NSGS_SHUFFLE_SEED] > 0) {
      srand((unsigned int)options->iparam[PLASTICITY_NSGS_SHUFFLE_SEED]);
    } else
      srand(1);
    scones = (unsigned int *)malloc(nc * sizeof(unsigned int));
    for (unsigned int i = 0; i < nc; ++i) {
      scones[i] = i;
    }
    uint_shuffle(scones, nc);
  }
  return scones;
}
static unsigned int *allocfreezingCones(PlasticityProblem *problem,
                                        SolverOptions *options) {
  unsigned int *fcones = 0;
  unsigned int nc = problem->numberOfCones;
  if (options->iparam[PLASTICITY_NSGS_FREEZING_CONE] > 0) {
    fcones = (unsigned int *)malloc(nc * sizeof(unsigned int));
    for (unsigned int i = 0; i < nc; ++i) {
      fcones[i] = 0;
    }
  }
  return fcones;
}

static int solveLocalReaction(UpdatePtr update_localproblem, SolverPtr local_solver,
                              CopyLocalReactionPtr copyLocalReaction, unsigned int cone,
                              PlasticityProblem *problem,
                              PlasticityProblem *localproblem, double *stress,
                              SolverOptions *localsolver_options, double localstress[3]) {
  (*update_localproblem)(cone, problem, localproblem, stress, localsolver_options);

  localsolver_options->iparam[PLASTICITY_CURRENT_CONE_NUMBER] = cone;

  copyLocalReaction(&(stress[cone * problem->dimension]), localstress);

  return (*local_solver)(localproblem, localstress, localsolver_options);
}

// static int file_exists(const char *fname) {
//   FILE *file;
//   if ((file = fopen(fname, "r"))) {
//     fclose(file);
//     return 1;
//   }
//   return 0;
// }

static void acceptLocalReactionFiltered(PlasticityProblem *localproblem,
                                        SolverOptions *localsolver_options, unsigned int cone,
                                        unsigned int iter, double *stress,
                                        double localstress[3]) {
  if (isnan(SOLVER_RESIDUAL(localsolver_options)) ||
      isinf(SOLVER_RESIDUAL(localsolver_options)) ||
      SOLVER_RESIDUAL(localsolver_options) > 1.0) {
    DEBUG_EXPR(plasticity2DProblem_display(localproblem));
    DEBUG_PRINTF(
        "Discard local stress for cone %i at iteration %i "
        "with local_error = %e\n",
        cone, iter, SOLVER_RESIDUAL(localsolver_options));

    /* #ifdef FCLIB_OUTPUT */

    /*     /\* printf("step counter value = %i\n", localsolver_options->iparam[19]); *\/ */
    /*     char fname[256]; */
    /*     fccounter++; */
    /*     snprintf(fname, sizeof(fname), "./local_problem/localproblem_%i_%i.hdf5", cone,
     */
    /*              localsolver_options->iparam[19]); */

    /*     if (file_exists(fname)) { */
    /*       /\* printf(" %s already dumped\n", fname); *\/ */
    /*     } else { */
    /*       printf("Dump %s\n", fname); */
    /*       int n = 100; */
    /*       char *title = (char *)malloc(n * sizeof(char)); */
    /*       strcpy(title, "Bad local problem dump in hdf5"); */
    /*       char *description = (char *)malloc(n * sizeof(char)); */
    /*       strcpy(description, "Rewriting in hdf5 from siconos "); */
    /*       strcat(description, fname); */
    /*       strcat(description, " in FCLIB format"); */
    /*       char *mathInfo = (char *)malloc(n * sizeof(char)); */
    /*       strcpy(mathInfo, "unknown"); */

    /*       plasticity2D_fclib_write(localproblem, title, description, mathInfo, fname, 3);
     */

    /*       printf("end of dump %s\n", fname); */
    /*       free(title); */
    /*       free(description); */
    /*       free(mathInfo); */
    /*     } */

    /* #endif */

    numerics_printf(
        "Discard local stress for cone %i at iteration %i "
        "with local_error = %e",
        cone, iter, SOLVER_RESIDUAL(localsolver_options));
  } else
    memcpy(&stress[cone * localproblem->dimension], localstress,
           sizeof(double) * localproblem->dimension);
}

static double calculateFullErrorAdaptiveInterval(PlasticityProblem *problem,
                                                 plasticity_2d_ComputeErrorPtr computeError,
                                                 SolverOptions *options, int iter,
                                                 double *stress, double *strainrate,
                                                 double tolerance, double norm_q) {
  double error = 1e+24;
  if (options->iparam[PLASTICITY_IPARAM_ERROR_EVALUATION_FREQUENCY] > 0) {
    if (iter % options->iparam[PLASTICITY_IPARAM_ERROR_EVALUATION_FREQUENCY] == 0) {
      (*computeError)(problem, stress, strainrate, tolerance, options, norm_q, &error);
      if (error > tolerance && options->iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] ==
                                   PLASTICITY_NSGS_ERROR_EVALUATION_ADAPTIVE)
        options->iparam[PLASTICITY_IPARAM_ERROR_EVALUATION_FREQUENCY] *= 2;
    }
    numerics_printf(
        "--------------- PLASTICITY_2D - NSGS - Iteration %i "
        "options->iparam[PLASTICITY_IPARAM_ERROR_EVALUATION_FREQUENCY] = %i, "
        "options->iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] = % i",
        iter, options->iparam[PLASTICITY_IPARAM_ERROR_EVALUATION_FREQUENCY],
        options->iparam[PLASTICITY_IPARAM_ERROR_EVALUATION]);
  } else
    (*computeError)(problem, stress, strainrate, tolerance, options, norm_q, &error);

  return error;
}

static double calculateFullErrorFinal(PlasticityProblem *problem, SolverOptions *options,
                                      plasticity_2d_ComputeErrorPtr computeError, double *stress,
                                      double *strainrate, double tolerance, double norm_q) {
  double absolute_error;
  (*computeError)(problem, stress, strainrate, tolerance, options, norm_q, &absolute_error);

  if (verbose > 0) {
    if (absolute_error > SOLVER_TOL(options)) {
      numerics_printf(
          "------- PLASTICITY_2D - NSGS - Warning absolute "
          "Residual = %14.7e is larger than required precision = %14.7e",
          absolute_error, SOLVER_TOL(options));
    } else {
      numerics_printf(
          "------- PLASTICITY_2D - NSGS - absolute "
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
    numerics_printf(
        "---- PLASTICITY_2D - NSGS - | %3d | %14.7e | %7.3e |",
        iter, error, tolerance);
  } else {
    numerics_printf(
        "---- PLASTICITY_2D - NSGS - | %3d | %14.7e | %7.3e |",
        iter, error, tolerance);
  }
  return hasNotConverged;
}

static int determine_convergence_with_full_final(PlasticityProblem *problem,
                                                 SolverOptions *options,
                                                 plasticity_2d_ComputeErrorPtr computeError,
                                                 double *stress, double *strainrate,
                                                 double *tolerance, double norm_q,
                                                 double error, int iter) {
  int hasNotConverged = 1;
  if (error < *tolerance) {
    hasNotConverged = 0;
    numerics_printf(
        "---- PLASTICITY_2D - NSGS - | %3d | %14.7e | %7.3e |",
        iter, error, *tolerance);

    double absolute_error =
        calculateFullErrorFinal(problem, options, computeError, stress, strainrate,
                                SOLVER_TOL(options), norm_q);
    if (absolute_error > SOLVER_TOL(options)) {
      if (error < DBL_EPSILON) {
        /* in this case, the relative error is very small
           (meaning that the nsgs loop does not
           improve accuracy).
           We try to tighten the local solver tolerance */
        SET_LOCAL_SOLVER_TOL(options->internalSolvers[0], fmax(LOCAL_SOLVER_TOL(options->internalSolvers[0]) / 100., DBL_EPSILON * 1e-6));
        numerics_printf(
            "------- PLASTICITY_2D - NSGS - We modify the local solver tolerance precision to reach "
            "accuracy to %e",
            LOCAL_SOLVER_TOL(options->internalSolvers[0]));

      } else {
        *tolerance = error / absolute_error * SOLVER_TOL(options);
	assert(*tolerance > 0.0 && "tolerance has to be positive");
        numerics_printf(
            "------- PLASTICITY_2D - NSGS - We modify the required incremental precision to reach "
            "accuracy to %e",
            *tolerance);
      }
      hasNotConverged = 1;
    } else {
      numerics_printf(
          "------- PLASTICITY_2D - NSGS - The incremental precision is sufficient to reach accuracy "
          "to %e",
          *tolerance);
    }

  } else {
    numerics_printf(
        "---- PLASTICITY_2D - NSGS - | %3d | %14.7e | %7.3e |",
        iter, error, *tolerance);
  }
  return hasNotConverged;
}

void plasticity_2d_nsgs(PlasticityProblem *problem, double *stress, double *strainrate, int *info,
               SolverOptions *options) {
  /* verbose=1; */


  /* plasticity_display(problem); */
  
  /* FILE *foutput = fopen("plasticity_2d_footing_100_theta0.05.dat", "w"); */
  /* // int info_output = */
  /* plasticity2D_printInFile(problem, foutput); */
  /* fclose(foutput); */

  /* int and double parameters */

  int *iparam = options->iparam;
  double *dparam = options->dparam;

  /* Number of cones */
  unsigned int nc = problem->numberOfCones;

  /* Maximum number of iterations */
  int itermax = SOLVER_MAX_ITER(options);

  /* Tolerance */
  double tolerance = SOLVER_TOL(options);
  SolverOptions *local_opts = options->internalSolvers[0];
  double local_tolerance_save = LOCAL_SOLVER_TOL(local_opts);
  double norm_q = cblas_dnrm2(nc * 3, problem->q, 1);
  double omega = dparam[PLASTICITY_NSGS_RELAXATION_VALUE];

  double norm_r[] = {1e24};
  if (options->numberOfInternalSolvers < 1) {
    *info = numerics_error("plasticity_2d_nsgs",
                   "The NSGS method needs options for the internal solvers, "
                   "options[0].numberOfInternalSolvers should be >= 1");
    return;
  }

  SolverOptions *localsolver_options = options->internalSolvers[0];
  plasticity_2d_ComputeErrorPtr computeError = NULL;

  struct LocalPLASTICITY_2DProblemFunctionToolkit *localProblemFunctionToolkit =
      localPLASTICITY_2DProblemFunctionToolkit_new();
  /* localPLASTICITY_2DProblemFunctionToolkit_display(localProblemFunctionToolkit); */

  PlasticityProblem *localproblem;

  double localstress[3];

  /*****  NSGS Iterations *****/
  int iter = 0;      /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  unsigned int cone; /* Number of the current row of blocks in M */
  unsigned int *scones = NULL;
  unsigned int *freeze_cones = NULL;

  /* Solver initialization continues below */

  SparseBlockStructuredMatrix *matrix1 = problem->M->matrix1;
  if (problem->M->storageType == NM_SPARSE) {
    if (problem->M->matrix1) {
      printf("Warning matrix 1 different from NULL");
    }

    problem->M->matrix1 = NM_extract_diagonal_blocks(problem->M, problem->dimension);
  }

  /*****  Initialize various solver options *****/
  localproblem = plasticity_2d_local_problem_allocate(problem);

  plasticity_2d_nsgs_initialize_local_solver(localProblemFunctionToolkit, &computeError, problem,
                                    localproblem, options);

  /* localProblemFunctionToolkit_display(localProblemFunctionToolkit); */
  scones = allocShuffledCones(problem, options);
  freeze_cones = allocfreezingCones(problem, options);
  /*****  Check solver options *****/
  if (!(iparam[PLASTICITY_NSGS_SHUFFLE] == PLASTICITY_NSGS_SHUFFLE_FALSE ||
        iparam[PLASTICITY_NSGS_SHUFFLE] == PLASTICITY_NSGS_SHUFFLE_TRUE ||
        iparam[PLASTICITY_NSGS_SHUFFLE] == PLASTICITY_NSGS_SHUFFLE_TRUE_EACH_LOOP)) {
    *info = numerics_error("plasticity_2d_nsgs",
                   "iparam[PLASTICITY_NSGS_SHUFFLE] must be equal to "
                   "PLASTICITY_NSGS_SHUFFLE_FALSE (0), "
                   "PLASTICITY_NSGS_SHUFFLE_TRUE (1) or "
                   "PLASTICITY_NSGS_SHUFFLE_TRUE_EACH_LOOP (2)");
    return;
  }

  if (!(iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] == PLASTICITY_NSGS_ERROR_EVALUATION_FULL ||
        iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] ==
            PLASTICITY_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL ||
        iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] == PLASTICITY_NSGS_ERROR_EVALUATION_LIGHT ||
        iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] ==
            PLASTICITY_NSGS_ERROR_EVALUATION_ADAPTIVE)) {
    *info = numerics_error("plasticity_2d_nsgs",
                   "iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] must be equal to "
                   "PLASTICITY_NSGS_ERROR_EVALUATION_FULL (0), "
                   "PLASTICITY_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL (1), "
                   "PLASTICITY_NSGS_ERROR_EVALUATION_LIGHT (2) or "
                   "PLASTICITY_NSGS_ERROR_EVALUATION_ADAPTIVE (3)");
    return;
  }

  /*****  NSGS Iterations *****/

  /* A special case for the most common options (should correspond
   * with mechanics_run.py **/
  if (iparam[PLASTICITY_NSGS_SHUFFLE] == PLASTICITY_NSGS_SHUFFLE_FALSE &&
      iparam[PLASTICITY_NSGS_FREEZING_CONE] == 0 &&
      iparam[PLASTICITY_NSGS_RELAXATION] == PLASTICITY_NSGS_RELAXATION_FALSE &&
      iparam[PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION] ==
          PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION_TRUE &&
      iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] == PLASTICITY_NSGS_ERROR_EVALUATION_LIGHT) {
    while ((iter < itermax) && (hasNotConverged > 0)) {
      ++iter;
      double light_error_sum = 0.0;

      plasticity_2d_set_internalsolver_tolerance(problem, options, localsolver_options, error);

      for (unsigned int i = 0; i < nc; ++i) {
        cone = i;

        solveLocalReaction(localProblemFunctionToolkit->update_local_problem,
                           localProblemFunctionToolkit->local_solver,
                           localProblemFunctionToolkit->copy_local_reaction, cone, problem,
                           localproblem, stress, localsolver_options, localstress);

        light_error_sum += localProblemFunctionToolkit->light_error_squared(
            localstress, &stress[cone * 3]);

        /* #if 0 */
        acceptLocalReactionFiltered(localproblem, localsolver_options, cone, iter, stress,
                                    localstress);
      }

      error = calculateLightError(light_error_sum, nc, stress, norm_r);

      hasNotConverged = determine_convergence(error, tolerance, iter, options);

      statsIterationCallback(problem, options, stress, strainrate, error);
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
      plasticity_2d_set_internalsolver_tolerance(problem, options, localsolver_options, error);

      unsigned int number_of_freezed_cone = 0;
      double tmp_criteria1 = tolerance * tolerance * 100 * 100;
      double tmp_criteria2 = *norm_r * *norm_r / (nc * nc * 1000);

      if (iparam[PLASTICITY_NSGS_FREEZING_CONE] > 0) {
        for (unsigned int i = 0; i < nc; ++i) {
          if (freeze_cones[i] > 0) number_of_freezed_cone++;
        }
        if (number_of_freezed_cone >= nc - 1) {
          // printf("number of freezed cone too large\n");
          for (unsigned int c = 0; c < nc; ++c) freeze_cones[c] = 0;
        }
      }
      for (unsigned int i = 0; i < nc; ++i) {
        if (iparam[PLASTICITY_NSGS_SHUFFLE] == PLASTICITY_NSGS_SHUFFLE_TRUE ||
            iparam[PLASTICITY_NSGS_SHUFFLE] == PLASTICITY_NSGS_SHUFFLE_TRUE_EACH_LOOP) {
          if (iparam[PLASTICITY_NSGS_SHUFFLE] == PLASTICITY_NSGS_SHUFFLE_TRUE_EACH_LOOP)
            uint_shuffle(scones, nc);
          cone = scones[i];
        } else
          cone = i;

        if (iparam[PLASTICITY_NSGS_FREEZING_CONE] > 0) {
          if (freeze_cones[cone] > 0) {
            /* we skip freeze cones */
            freeze_cones[cone] -= 1;
            continue;
          }
        }

        solveLocalReaction(localProblemFunctionToolkit->update_local_problem,
                           localProblemFunctionToolkit->local_solver,
                           localProblemFunctionToolkit->copy_local_reaction, cone, problem,
                           localproblem, stress, localsolver_options, localstress);

        if (iparam[PLASTICITY_NSGS_RELAXATION] == PLASTICITY_NSGS_RELAXATION_TRUE)
          localProblemFunctionToolkit->perform_relaxation(localstress, &stress[cone * 3],
                                                          omega);

        light_error_2 = localProblemFunctionToolkit->light_error_squared(localstress,
                                                                         &stress[cone * 3]);

        light_error_sum += light_error_2;

        /* int test =100; */
        /* if (cone == test) */
        /* { */
        /*   printf("stress[%i] = %16.8e\t",3*cone-1,stress[3*cone]); */
        /*   printf("localstress[%i] = %16.8e\n",2,localstress[0]); */
        /* } */

        if (iparam[PLASTICITY_NSGS_FREEZING_CONE] > 0) {
          double squared_norm_localstress =
              localProblemFunctionToolkit->squared_norm(localstress);
          int relative_convergence_criteria =
              light_error_2 <= tmp_criteria1 * squared_norm_localstress;
          int small_stress_criteria = squared_norm_localstress <= tmp_criteria2;

          if ((relative_convergence_criteria || small_stress_criteria) && iter >= 10)
          /* if ((light_error_2 *squared_norm(localstress) <= tolerance*tolerance/(nc*nc*10)
           */
          /*      || squared_norm(localstress) <=  (*norm_r* *norm_r/(nc*nc*1000))) */
          /*     && iter >=10) */
          {
            /* we  freeze the cone for n iterations*/
            freeze_cones[cone] = iparam[PLASTICITY_NSGS_FREEZING_CONE];

            DEBUG_EXPR(printf("first criteria : light_error_2*squared_norm(localstress) <= "
                              "tolerance*tolerance/(nc*nc*10) ==> %e <= %e, bool =%i\n",
                              light_error_2 * squared_norm(localstress),
                              tolerance * tolerance / (nc * nc * 10),
                              relative_convergence_criteria);
                       printf("second criteria :  squared_norm(localstress) <=  (*norm_r* "
                              "*norm_r/(nc*nc))/1000. ==> %e <= %e, bool =%i \n",
                              squared_norm(localstress),
                              *norm_r * *norm_r / (nc * nc * 1000), small_stress_criteria);
                       printf("Cone % i is freezed for %i steps\n", cone,
                              iparam[PLASTICITY_NSGS_FREEZING_CONE]););
          }
        }

        if (iparam[PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION] ==
            PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION_TRUE)
          acceptLocalReactionFiltered(localproblem, localsolver_options, cone, iter, stress,
                                      localstress);
        else
          acceptLocalReactionUnconditionally(cone, stress, localstress);
      }

      /* DEBUG_EXPR( */
      /*   if(iparam[PLASTICITY_NSGS_FREEZING_CONE] >0) */
      /*   { */
      /*     int frozen_cone=0; */
      /*     for(unsigned int ii = 0 ; ii < nc ; ++ii) if (freeze_cones[ii] >0)
       * frozen_cone++; */
      /*     numerics_printf_verbose(1,"number of frozen cones %i at iter : %i",
       * frozen_cone, iter ); */
      /*   } */
      /*   ); */

      if (iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] ==
          PLASTICITY_NSGS_ERROR_EVALUATION_LIGHT) {
        error = calculateLightError(light_error_sum, nc, stress, norm_r);
        hasNotConverged = determine_convergence(error, tolerance, iter, options);
      } else if (iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] ==
                 PLASTICITY_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
        error = calculateLightError(light_error_sum, nc, stress, norm_r);
        hasNotConverged =
            determine_convergence_with_full_final(problem, options, computeError, stress,
                                                  strainrate, &tolerance, norm_q, error, iter);

        if (!(tolerance > 0.0)) {
          numerics_warning("plasticity_2d_nsgs", "tolerance has to be positive!!");
          numerics_warning("plasticity_2d_nsgs", "we stop the iterations");
          break;
        }

      } else if (iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] ==
                 PLASTICITY_NSGS_ERROR_EVALUATION_FULL) {
        error = calculateFullErrorAdaptiveInterval(problem, computeError, options, iter,
                                                   stress, strainrate, tolerance, norm_q);
        hasNotConverged = determine_convergence(error, tolerance, iter, options);
      }

      statsIterationCallback(problem, options, stress, strainrate, error);

      /* if(iparam[PLASTICITY_NSGS_FREEZING_CONE] >0) */
      /* { */
      /*   int frozen_cone=0; */
      /*   for(unsigned int i = 0 ; i < nc ; ++i) */
      /*   { */
      /*     if (freeze_cones[i] >0) */
      /*     { */
      /*       frozen_cone++; */
      /*     } */
      /*   } */
      /*   printf("number of frozen cones %i at iter : %i over number of cones: %i\n",
       * frozen_cone, iter, nc ); */
      /* } */
    }
  }

  /* Full criterium */
  if (iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] ==
      PLASTICITY_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
    error = calculateFullErrorFinal(problem, options, computeError, stress, strainrate,
                                    tolerance, norm_q);

    hasNotConverged = determine_convergence(error, tolerance, iter, options);
  }

  *info = hasNotConverged;

  /** return parameter values */
  /* dparam[SICONOS_DPARAM_TOL] = tolerance; */
  SET_SOLVER_RESIDUAL(options, error);
  SET_SOLVER_ITER_DONE(options, iter);
  SET_LOCAL_SOLVER_TOL(local_opts, local_tolerance_save);

  /** Free memory **/

  if (problem->M->storageType == NM_SPARSE) {
    SBM_clear_block(problem->M->matrix1);
    SBM_clear(problem->M->matrix1);
    problem->M->matrix1 = matrix1;
  }
  localProblemFunctionToolkit->free_local_solver(problem, localproblem, localsolver_options);

  plasticity_2d_local_problem_free(localproblem, problem);
  if (scones) free(scones);
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers PLASTICITY_2D_NSGS in the global solver registry, enabling:
 * - Dynamic solver lookup by ID
 * - Runtime solver introspection
 * - Elimination of giant switch statements in drivers
 */

static int plasticity_2d_nsgs_init_wrap(void* problem, SolverOptions* options) {
  plasticity_2d_nsgs_set_default(options);
  return NUMERICS_OK;
}

static int plasticity_2d_nsgs_solve_wrap(void* problem, double* reaction,
                                double* velocity, SolverOptions* options) {
  int info = NUMERICS_OK;
  plasticity_2d_nsgs((PlasticityProblem*)problem, reaction, velocity, &info, options);
  return info;
}

static void plasticity_2d_nsgs_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(PLASTICITY_2D_NSGS, "PLASTICITY_2D_NSGS",
                "Non-smooth Gauss-Seidel for 2D Mohr Coulomb Plasticity",
                plasticity_2d_nsgs_init_wrap,
                plasticity_2d_nsgs_solve_wrap,
                plasticity_2d_nsgs_free_wrap,
                NULL,  /* error function */
                plasticity_2d_nsgs_set_default,  /* set_default */
                1000,  /* default_max_iter */
                1e-4,  /* default_tol */
                0      /* is_local_solver */);

void plasticity_2d_nsgs_set_default(SolverOptions *options) {
  options->iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] =
      PLASTICITY_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL;
  // options->iparam[PLASTICITY_IPARAM_ERROR_EVALUATION] =
  // PLASTICITY_NSGS_ERROR_EVALUATION_FULL;
  options->iparam[PLASTICITY_IPARAM_INTERNAL_ERROR_STRATEGY] =
      PLASTICITY_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE;
  /* options->iparam[PLASTICITY_IPARAM_INTERNAL_ERROR_STRATEGY] =
   * PLASTICITY_INTERNAL_ERROR_STRATEGY_ADAPTIVE; */
  /* options->iparam[PLASTICITY_IPARAM_INTERNAL_ERROR_STRATEGY] =
   * PLASTICITY_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONE; */
  options->iparam[PLASTICITY_NSGS_SHUFFLE] = PLASTICITY_NSGS_SHUFFLE_FALSE;
  options->iparam[PLASTICITY_NSGS_SHUFFLE_SEED] = 0;
  options->iparam[PLASTICITY_NSGS_FREEZING_CONE] = 0;
  options->iparam[PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION] =
      PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION_FALSE;
  options->iparam[PLASTICITY_NSGS_RELAXATION] = PLASTICITY_NSGS_RELAXATION_FALSE;
  options->iparam[PLASTICITY_IPARAM_ERROR_EVALUATION_FREQUENCY] = 0;
  SOLVER_TOL(options) = 1e-4;
  options->dparam[PLASTICITY_DPARAM_INTERNAL_ERROR_RATIO] = 10.0;
  // Internal solver - allocate if needed
  if (options->numberOfInternalSolvers == 0) {
    options->numberOfInternalSolvers = 1;
    options->internalSolvers = calloc(1, sizeof(SolverOptions*));
  }
  assert(options->numberOfInternalSolvers == 1);
  options->internalSolvers[0] = solver_options_create(PLASTICITY_2D_ONECONE_NSN_GP_HYBRID);
}
