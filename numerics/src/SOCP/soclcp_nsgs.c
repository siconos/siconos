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
#include <math.h>    // for pow, sqrt
#ifndef __cplusplus
#include <stdbool.h>  // for false
#endif
#include <stdio.h>   // for printf, NULL
#include <stdlib.h>  // for malloc, free
#include <string.h>  // for memcpy

#include "NSSTools.h"                                     // for max
#include "NumericsFwd.h"                                  // for SecondOrder...
#include "NumericsMatrix.h"                               // for NumericsMatrix
#include "SOCLCP_Solvers.h"                               // for FreeSolverN...
#include "SOCLCP_cst.h"                                   // for SICONOS_SOC...
#include "SecondOrderConeLinearComplementarityProblem.h"  // for SecondOrder...
#include "SolverOptions.h"                                // for SolverOptions
#include "numerics_verbose.h"                             // for verbose
#include "soclcp_compute_error.h"                         // for soclcp_comp...
#include "soclcp_projection.h"                            // for soclcp_proj...
#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void soclcp_nsgs_update(int cone, SecondOrderConeLinearComplementarityProblem* problem,
                        SecondOrderConeLinearComplementarityProblem* localproblem, double* r,
                        SolverOptions* options) {
  /* Build a local problem for a specific cone
     r corresponds to the global vector (size n) of the global problem.
  */
  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific cone
   r corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  soclcp_nsgs_fillMLocal(problem, localproblem, cone);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products
     MLocal.rBlock, excluding the block corresponding to the current cone. ****/
  soclcp_nsgs_computeqLocal(problem, localproblem, r, cone, options);

  /* coefficient for current block*/
  localproblem->tau[0] = problem->tau[cone];

  /* index for current block*/
  localproblem->coneIndex[0] = 0;

  /* coefficient for current block*/
  localproblem->n = problem->coneIndex[cone + 1] - problem->coneIndex[cone];
}
void soclcp_initializeLocalSolver_nsgs(
    Solver_soclcp_Ptr* solve, Update_soclcp_Ptr* update, FreeSolverNSGS_soclcp_Ptr* freeSolver,
    ComputeError_soclcp_Ptr* computeError,
    SecondOrderConeLinearComplementarityProblem* problem,
    SecondOrderConeLinearComplementarityProblem* localproblem,
    SolverOptions* localsolver_options) {
  /** Connect to local solver */
  switch (localsolver_options->solverId) {
    /* Projection */
    case SICONOS_SOCLCP_ProjectionOnCone: {
      *solve = &soclcp_projectionOnCone_solve;
      *update = &soclcp_nsgs_update;
      *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_projection_free;
      *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error;
      soclcp_projection_initialize(problem, localproblem, localsolver_options);
      break;
    }
    case SICONOS_SOCLCP_ProjectionOnConeWithLocalIteration: {
      *solve = &soclcp_projectionOnConeWithLocalIteration_solve;
      *update = &soclcp_nsgs_update;
      *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_projectionOnConeWithLocalIteration_free;
      *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error;
      soclcp_projectionOnConeWithLocalIteration_initialize(problem, localproblem,
                                                           localsolver_options);
      break;
    }
    case SICONOS_SOCLCP_ProjectionOnConeWithRegularization: {
      *solve = &soclcp_projectionOnCone_solve;
      *update = &soclcp_projection_update_with_regularization;
      *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_projection_with_regularization_free;
      *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error;
      soclcp_projection_initialize_with_regularization(problem, localproblem);
      break;
    }
    /* /\* Newton solver (Alart-Curnier) *\/ */
    /* case SICONOS_SOCLCP_AlartCurnierNewton: */
    /* { */
    /*   *solve = &soclcp_Newton_solve; */
    /*   *update = &soclcp_AC_update; */
    /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_Newton_free; */
    /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
    /*   soclcp_Newton_initialize(problem, localproblem, localsolver_options); */
    /*   break; */
    /* } */
    /* case SICONOS_SOCLCP_DampedAlartCurnierNewton: */
    /* { */
    /*   *solve = &soclcp_Newton_solve; */
    /*   *update = &soclcp_AC_update; */
    /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_Newton_free; */
    /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
    /*   soclcp_Newton_initialize(problem, localproblem, localsolver_options); */
    /*   break; */
    /* } */
    /* /\* Newton solver (Glocker-Fischer-Burmeister)*\/ */
    /* case SICONOS_SOCLCP_NCPGlockerFBNewton: */
    /* { */
    /*   *solve = &soclcp_Newton_solve; */
    /*   *update = &NCPGlocker_update; */
    /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_Newton_free; */
    /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
    /*   // *computeError = &fake_compute_error; */
    /*   soclcp_Newton_initialize(problem, localproblem, localsolver_options); */
    /*   break; */
    /* } */
    /* /\* Path solver (Glocker Formulation) *\/ */
    /* case SICONOS_SOCLCP_NCPGlockerFBPATH: */
    /* { */
    /*   *solve = &soclcp_Path_solve; */
    /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_Path_free; */
    /*   *update = &NCPGlocker_update; */
    /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
    /*   // *computeError = &fake_compute_error; */
    /*   soclcp_Path_initialize(problem, localproblem, localsolver_options); */
    /*   break; */
    /* } */

    /* /\* Fixed Point solver (Glocker Formulation) *\/ */
    /* case SICONOS_SOCLCP_NCPGlockerFBFixedPoint: */
    /* { */
    /*   *solve = &soclcp_FixedP_solve; */
    /*   *update = &NCPGlocker_update; */
    /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_FixedP_free; */
    /*   *computeError = &fake_compute_error_nsgs; */
    /*   soclcp_FixedP_initialize(problem, localproblem, localsolver_options); */
    /*   break; */
    /* } */
    /* /\*iparam[4] > 10 are reserved for Tresca resolution *\/ */
    /* case SICONOS_SOCLCP_projectionOnCylinder: */
    /* { */
    /*   *solve = &soclcp_projectionOnCylinder_solve; */
    /*   *update = &soclcp_projectionOnCylinder_update; */
    /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_projection_free; */
    /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_Tresca_compute_error; */
    /*   soclcp_projection_initialize(problem, localproblem); */
    /*   break; */
    /* } */
    /* case SICONOS_SOCLCP_QUARTIC: */
    /* { */
    /*   *solve = &soclcp_unitary_enumerative_solve; */
    /*   *update = &soclcp_nsgs_update; */
    /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_unitary_enumerative_free; */
    /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
    /*   soclcp_unitary_enumerative_initialize(localproblem); */
    /*   break; */
    /* } */
    /* case SICONOS_SOCLCP_QUARTIC_NU: */
    /* { */
    /*   *solve = &soclcp_unitary_enumerative_solve; */
    /*   *update = &soclcp_nsgs_update; */
    /*   *freeSolver = (FreeSolverNSGS_soclcp_Ptr)&soclcp_unitary_enumerative_free; */
    /*   *computeError = (ComputeError_soclcp_Ptr)&soclcp_compute_error; */
    /*   soclcp_unitary_enumerative_initialize(localproblem); */
    /*   break; */
    /* } */
    default: {
      fprintf(stderr, "Numerics, soclcp_nsgs failed. Unknown internal solver : %s.\n",
              solver_options_id_to_name(localsolver_options->solverId));
      exit(EXIT_FAILURE);
    }
  }
}
void soclcp_nsgs_computeqLocal(SecondOrderConeLinearComplementarityProblem* problem,
                               SecondOrderConeLinearComplementarityProblem* localproblem,
                               double* r, int cone, SolverOptions* options) {
  double* qLocal = localproblem->q;
  int n = problem->n;
  int normal = problem->coneIndex[cone];
  int dim = problem->coneIndex[cone + 1] - problem->coneIndex[cone];

  /* r current block set to zero, to exclude current cone block */
  int i;
  double* rsave = options->dWork;

  /* qLocal computation*/
  for (i = 0; i < dim; i++) {
    qLocal[i] = problem->q[normal + i];
  }

  NM_row_prod_no_diag(n, dim, cone, normal, problem->M, r, qLocal, rsave, false);
}

void soclcp_nsgs_fillMLocal(SecondOrderConeLinearComplementarityProblem* problem,
                            SecondOrderConeLinearComplementarityProblem* localproblem,
                            int cone) {
  int coneStart = problem->coneIndex[cone];
  NM_extract_diag_block(problem->M, cone, coneStart, problem->coneIndex[cone + 1] - coneStart,
                        &localproblem->M->matrix0);
}

/* swap two indices */
void uint_swap(unsigned int* a, unsigned int* b);

/* shuffle an unsigned array */
void uint_shuffle(unsigned int* a, unsigned int n);

void soclcp_nsgs(SecondOrderConeLinearComplementarityProblem* problem, double* r, double* v,
                 int* info, SolverOptions* options) {
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of cones */
  unsigned int nc = problem->nc;
  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];

  if (*info == 0) return;

  if (options->numberOfInternalSolvers < 1) {
    numerics_error("soclcp_nsgs",
                   "The NSGS method needs options for the internal solvers, "
                   "options[0].numberOfInternalSolvers should be >1");
  }
  assert(&options[1]);

  SolverOptions* localsolver_options = options->internalSolvers[0];

  Solver_soclcp_Ptr local_solver = NULL;
  Update_soclcp_Ptr update_localproblem = NULL;
  FreeSolverNSGS_soclcp_Ptr freeSolver = NULL;
  ComputeError_soclcp_Ptr computeError = NULL;

  unsigned int cone;
  int isConeDimensionsEqual = 1;
  /* Connect local solver and local problem*/
  unsigned int dim_max = problem->coneIndex[1] - problem->coneIndex[0];
  for (cone = 1; cone < nc; cone++) {
    if (dim_max != problem->coneIndex[cone + 1] - problem->coneIndex[cone]) {
      isConeDimensionsEqual = 0;
      dim_max = max(dim_max, problem->coneIndex[cone + 1] - problem->coneIndex[cone]);
    }
  }

  SecondOrderConeLinearComplementarityProblem* localproblem;

  if (isConeDimensionsEqual) {
    localproblem = (SecondOrderConeLinearComplementarityProblem*)malloc(
        sizeof(SecondOrderConeLinearComplementarityProblem));
    localproblem->nc = 1;
    localproblem->n = dim_max;
    localproblem->q = (double*)malloc(dim_max * sizeof(double));
    localproblem->tau = (double*)malloc(sizeof(double));
    localproblem->coneIndex = (unsigned int*)malloc(2 * sizeof(unsigned int));
    localproblem->coneIndex[0] = 0;
    localproblem->coneIndex[1] = dim_max;

    if (problem->M->storageType != NM_SPARSE_BLOCK) {
      localproblem->M = NM_create_from_data(NM_DENSE, dim_max, dim_max,
                                            malloc(dim_max * dim_max * sizeof(double)));
    } else {
      localproblem->M = NM_new();
    }
  } else {
    fprintf(stderr, "soclcp_nsgs error: not yet implemented.\n");
    exit(EXIT_FAILURE);
  }

  soclcp_initializeLocalSolver_nsgs(&local_solver, &update_localproblem,
                                    (FreeSolverNSGS_soclcp_Ptr*)&freeSolver, &computeError,
                                    problem, localproblem, localsolver_options);

  /*****  NSGS Iterations *****/
  int iter = 0;      /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  unsigned int* scones = NULL;

  if (iparam[SICONOS_IPARAM_NSGS_SHUFFLE]) /* shuffle */
  {
    scones = (unsigned int*)malloc(nc * sizeof(unsigned int));
    for (unsigned int i = 0; i < nc; ++i) {
      scones[i] = i;
    }
    uint_shuffle(scones, nc);
  }

  if (iparam[SICONOS_IPARAM_ERROR_EVALUATION] == SICONOS_ERROR_LIGHT_EVALUATION ||
      iparam[SICONOS_IPARAM_ERROR_EVALUATION] == SICONOS_ERROR_LIGHT_EVALUATION_NO_UPDATE) {
    int n = problem->n;
    double* rold =
        (double*)malloc(n * sizeof(double));  // save memory if isConeDimensionsEqual
    while ((iter < itermax) && (hasNotConverged > 0)) {
      ++iter;
      /* Loop through the cone  */
      // cblas_dcopy( n , q , incx , v , incy );
      error = 0.0;
      memcpy(rold, r, n * sizeof(double));
      for (unsigned int i = 0; i < nc; ++i) {
        if (iparam[SICONOS_IPARAM_NSGS_SHUFFLE]) {
          cone = scones[i];
        } else {
          cone = i;
        }

        if (verbose > 1) printf("--------------- Cone Number %i\n", cone);
        (*update_localproblem)(cone, problem, localproblem, r, localsolver_options);

        localsolver_options->iparam[SICONOS_IPARAM_SOCLCP_PROJECTION_CONE_INDEX] = cone;

        (*local_solver)(localproblem, &(r[problem->coneIndex[cone]]), localsolver_options);
      }
      for (int i = 0; i < n; i++) {
        error += pow(r[i] - rold[i], 2);
      }

      /* **** Criterium convergence **** */
      error = sqrt(error);
      if (verbose > 0)
        printf("--------------- SOCLP - NSGS - Iteration %i Residual = %14.7e >= %7.4e\n",
               iter, error, options->dparam[SICONOS_DPARAM_TOL]);
      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }

    if (iparam[SICONOS_IPARAM_ERROR_EVALUATION] ==
        SICONOS_ERROR_LIGHT_EVALUATION) /* Full criterium */
    {
      double absolute_error;
      (*computeError)(problem, r, v, tolerance, options, &absolute_error);
      if (verbose > 0) {
        if (absolute_error > error) {
          printf(
              "--------------- SOCLCP - NSGS - Warning absolute Residual = %14.7e is larger "
              "than incremental error = %14.7e\n",
              absolute_error, error);
        }
      }
    }
    free(rold);
  } else {
    if (iparam[SICONOS_IPARAM_NSGS_SHUFFLE]) {
      int withRelaxation = iparam[SICONOS_IPARAM_SOCLCP_NSGS_WITH_RELAXATION];
      if (withRelaxation) {
        int n = problem->n;
        double* rold =
            (double*)malloc(n * sizeof(double));  // save memory if isConeDimensionsEqual
        double omega = dparam[SICONOS_DPARAM_SOCLCP_NSGS_RELAXATION];
        unsigned int dim;
        while ((iter < itermax) && (hasNotConverged > 0)) {
          ++iter;
          memcpy(rold, r, n * sizeof(double));
          /* Loop through the cone points */
          // cblas_dcopy( n , q , incx , v , incy );
          for (unsigned int i = 0; i < nc; ++i) {
            cone = scones[i];

            if (verbose > 1) printf("--------------- Cone Number %i\n", cone);
            (*update_localproblem)(cone, problem, localproblem, r, localsolver_options);
            localsolver_options->iparam[SICONOS_IPARAM_SOCLCP_PROJECTION_CONE_INDEX] = cone;
            (*local_solver)(localproblem, &(r[problem->coneIndex[cone]]), localsolver_options);

            dim = problem->coneIndex[cone + 1] - problem->coneIndex[cone];
            for (unsigned int i = 0; i < dim; ++i) {
              r[problem->coneIndex[cone] + i] =
                  omega * r[problem->coneIndex[cone] + i] +
                  (1.0 - omega) * rold[problem->coneIndex[cone] + i];
            }
          }

          /* **** Criterium convergence **** */
          (*computeError)(problem, r, v, tolerance, options, &error);

          if (verbose > 0)
            printf("--------------- SOCLP - NSGS - Iteration %i Residual = %14.7e >= %7.4e\n",
                   iter, error, options->dparam[SICONOS_DPARAM_TOL]);

          if (error < tolerance) hasNotConverged = 0;
          *info = hasNotConverged;

          if (options->callback) {
            options->callback->collectStatsIteration(options->callback->env, problem->n, r, v,
                                                     error, NULL);
          }
        }
        free(rold);
      } else {
        while ((iter < itermax) && (hasNotConverged > 0)) {
          ++iter;
          /* Loop through the cone  */
          // cblas_dcopy( n , q , incx , v , incy );
          for (unsigned int i = 0; i < nc; ++i) {
            cone = scones[i];

            if (verbose > 1) printf("--------------- Cone Number %i\n", cone);
            (*update_localproblem)(cone, problem, localproblem, r, localsolver_options);
            localsolver_options->iparam[SICONOS_IPARAM_SOCLCP_PROJECTION_CONE_INDEX] = cone;
            (*local_solver)(localproblem, &(r[problem->coneIndex[cone]]), localsolver_options);
          }

          /* **** Criterium convergence **** */
          (*computeError)(problem, r, v, tolerance, options, &error);

          if (verbose > 0)
            printf("--------------- SOCLP - NSGS - Iteration %i Residual = %14.7e >= %7.4e\n",
                   iter, error, options->dparam[SICONOS_DPARAM_TOL]);

          if (error < tolerance) hasNotConverged = 0;
          *info = hasNotConverged;

          if (options->callback) {
            options->callback->collectStatsIteration(options->callback->env, problem->n, r, v,
                                                     error, NULL);
          }
        }
      }

    } else {
      int withRelaxation = iparam[SICONOS_IPARAM_SOCLCP_NSGS_WITH_RELAXATION];
      if (withRelaxation) {
        int n = problem->n;
        unsigned int dim;
        double* rold =
            (double*)malloc(n * sizeof(double));  // save memory if isConeDimensionsEqual
        double omega = dparam[SICONOS_DPARAM_SOCLCP_NSGS_RELAXATION];
        while ((iter < itermax) && (hasNotConverged > 0)) {
          ++iter;
          for (int i = 0; i < n; i++) rold[i] = r[i];

          /* Loop through the cones  */
          // cblas_dcopy( n , q , incx , v , incy );
          for (cone = 0; cone < nc; ++cone) {
            if (verbose > 1) printf("--------------- Cone Number %i\n", cone);
            (*update_localproblem)(cone, problem, localproblem, r, localsolver_options);
            localsolver_options->iparam[SICONOS_IPARAM_SOCLCP_PROJECTION_CONE_INDEX] = cone;
            (*local_solver)(localproblem, &(r[problem->coneIndex[cone]]), localsolver_options);

            dim = problem->coneIndex[cone + 1] - problem->coneIndex[cone];
            for (unsigned int i = 0; i < dim; ++i) {
              r[problem->coneIndex[cone] + i] =
                  omega * r[problem->coneIndex[cone] + i] +
                  (1.0 - omega) * rold[problem->coneIndex[cone] + i];
            }
          }

          /* **** Criterium convergence **** */
          (*computeError)(problem, r, v, tolerance, options, &error);

          if (verbose > 0)
            printf("--------------- SOCLP - NSGS - Iteration %i Residual = %14.7e >= %7.4e\n",
                   iter, error, options->dparam[SICONOS_DPARAM_TOL]);

          if (error < tolerance) hasNotConverged = 0;
          *info = hasNotConverged;

          if (options->callback) {
            options->callback->collectStatsIteration(options->callback->env, problem->n, r, v,
                                                     error, NULL);
          }
        }
        free(rold);

      } else {
        while ((iter < itermax) && (hasNotConverged > 0)) {
          ++iter;
          /* Loop through the cones  */
          // cblas_dcopy( n , q , incx , v , incy );
          for (cone = 0; cone < nc; ++cone) {
            if (verbose > 1) printf("--------------- Cone Number %i\n", cone);
            (*update_localproblem)(cone, problem, localproblem, r, localsolver_options);
            localsolver_options->iparam[SICONOS_IPARAM_SOCLCP_PROJECTION_CONE_INDEX] = cone;
            (*local_solver)(localproblem, &(r[problem->coneIndex[cone]]), localsolver_options);
          }

          /* **** Criterium convergence **** */
          (*computeError)(problem, r, v, tolerance, options, &error);

          if (verbose > 0)
            printf("--------------- SOCLP - NSGS - Iteration %i Residual = %14.7e >= %7.4e\n",
                   iter, error, options->dparam[SICONOS_DPARAM_TOL]);

          if (error < tolerance) hasNotConverged = 0;
          *info = hasNotConverged;

          if (options->callback) {
            options->callback->collectStatsIteration(options->callback->env, problem->n, r, v,
                                                     error, NULL);
          }
        }
      }
    }
  }
  dparam[SICONOS_DPARAM_TOL] = tolerance;
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

  /***** Free memory *****/
  (*freeSolver)(problem, localproblem, localsolver_options);

  if (problem->M->storageType == NM_SPARSE_BLOCK) localproblem->M->matrix0 = NULL;

  freeSecondOrderConeLinearComplementarityProblem(localproblem);

  if (scones) /* shuffle */
  {
    free(scones);
  }
}

void soclcp_nsgs_set_default(SolverOptions* options) {
  options->iparam[SICONOS_IPARAM_ERROR_EVALUATION] = SICONOS_ERROR_FULL_EVALUATION;
  options->iparam[SICONOS_IPARAM_SOCLCP_NSGS_WITH_RELAXATION] = 0;
  options->iparam[SICONOS_IPARAM_NSGS_SHUFFLE] = 0;
  options->dparam[SICONOS_DPARAM_SOCLCP_NSGS_RELAXATION] = 1.;
  assert(options->numberOfInternalSolvers == 1);
  options->internalSolvers[0] =
      solver_options_create(SICONOS_SOCLCP_ProjectionOnConeWithLocalIteration);
  options->internalSolvers[0]->dparam[SICONOS_DPARAM_SOCLCP_PROJECTION_RHO] = 0.;
}
