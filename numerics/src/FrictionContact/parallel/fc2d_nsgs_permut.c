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
#include <math.h>    // for fabs, sqrt, INFINITY
#include <stdio.h>   // for NULL, fprintf, printf, stderr
#include <stdlib.h>  // for free, malloc, calloc

#include "FrictionContactProblem.h"        // for FrictionContactProblem
#include "FrictionContact_options.h"       // for SICONOS_FRICTION_3D_IPARAM...
#include "LCP_Solvers.h"                   // for lcp_nsgs_SBM_buildLocalPro...
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NumericsFwd.h"                   // for SolverOptions, LinearCompl...
#include "NumericsMatrix.h"                // for NumericsMatrix, RawNumeric...
#include "SiconosBlas.h"                   // for cblas_dnrm2
#include "SolverOptions.h"                 // for SolverOptions, SICONOS_DPA...
#include "SparseBlockMatrix.h"             // for SparseBlockStructuredMatrix
#include "fc2d_Solvers.h"                  // for fc2d_nsgs_sbm, fc2d_spa...
#include "fc2d_compute_error.h"            // for fc2d_compute_error
#include "graph_tools.h"
#include "numerics_verbose.h"  // for numerics_printf, verbose

/* Solver registration system */
#include "numerics_errors.h"
#include "solver_registry.h"

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES 1 */
#include "siconos_debug.h"  // for DEBUG_BEGIN, DEBUG_END
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

#define SGN(x) ((x) < 0 ? -1 : (x) > 0 ? 1 : 0)

static SparseBlockCoordinateMatrix* fc3d_extract_diagonal_blocks(
    FrictionContactProblem* problem) {
  unsigned int nc = (unsigned int)problem->numberOfContacts;

  SparseBlockCoordinateMatrix* sbcm = SBCM_new();
  sbcm->blocknumber0 = nc;
  sbcm->blocknumber1 = nc;
  unsigned int nbblocks = nc;
  sbcm->nbblocks = nc;
  /* sbcm->row = (unsigned int *)  malloc(sizeof(unsigned int) *   nbblocks); */
  /* sbcm->column = (unsigned int *)  malloc(sizeof(unsigned int) *  nbblocks); */
  /* for(unsigned int i = 0; i < sbcm->nbblocks; ++i) */
  /* { */
  /*   sbcm->row[i] = i; */
  /* } */
  /* for(unsigned int i = 0; i < sbcm->nbblocks; ++i) */
  /* { */
  /*   sbcm->column[i] = i; */
  /* } */
  sbcm->block = (double**)malloc(sizeof(double*) * nbblocks);
  //  sbcm->blocksize0 = (unsigned int *) malloc(sizeof(unsigned int) * nbblocks);
  //  sbcm->blocksize1 = (unsigned int *) malloc(sizeof(unsigned int) * nbblocks);

  for (unsigned int contact = 0; contact < nc; ++contact) {
    //    sbcm->blocksize0[contact]=2;
    //    sbcm->blocksize1[contact]=2;
    if (problem->M->storageType != NM_SPARSE_BLOCK)
      sbcm->block[contact] = (double*)calloc(sizeof(double*), 4);
    NM_extract_diag_block2(problem->M, (int)contact, &sbcm->block[contact]);
  }

  return sbcm;
}
static SparseBlockCoordinateMatrix* fc3d_free_diagonal_blocks(
    FrictionContactProblem* problem, SparseBlockCoordinateMatrix* sbcm) {
  if (sbcm->row) free(sbcm->row);
  if (sbcm->column) free(sbcm->column);
  for (unsigned int contact = 0; contact < sbcm->nbblocks; ++contact) {
    if (problem->M->storageType != NM_SPARSE_BLOCK) {
      if (sbcm->block[contact]) free(sbcm->block[contact]);
    }
  }
  if (sbcm->blocksize0) free(sbcm->blocksize0);
  if (sbcm->blocksize0) free(sbcm->blocksize0);
  if (sbcm->block) free(sbcm->block);
  free(sbcm);
  return NULL;
}

static void fc2d_nsgs_buildLocalProblem(int contact, FrictionContactProblem* problem,
                                        SparseBlockCoordinateMatrix* diagonal_blocks,
                                        LinearComplementarityProblem* local_problem,
                                        double* reaction) {
  // NM_extract_diag_block2(problem->M, contact, &local_problem->M->matrix0);
  local_problem->M->matrix0 = diagonal_blocks->block[contact];

  local_problem->M->size0 = 2;  // Necessary ?
  local_problem->M->size1 = 2;

  local_problem->q[0] = problem->q[contact * 2];
  local_problem->q[1] = problem->q[contact * 2 + 1];
  NM_row_prod_no_diag2(2 * (size_t)problem->numberOfContacts, contact, 2 * (size_t)contact,
                       problem->M, reaction, local_problem->q, false);

  DEBUG_EXPR(NM_display(local_problem->M););
  DEBUG_EXPR(NV_display(local_problem->q, 2););
}

static inline double light_error_squared(double localreaction[2], double* oldreaction) {
  double x0 = oldreaction[0] - localreaction[0];
  double x1 = oldreaction[1] - localreaction[1];
  return x0 * x0 + x1 * x1;
}
static inline double squared_norm(double localreaction[2]) {
  return (localreaction[0] * localreaction[0] + localreaction[1] * localreaction[1]);
}

static inline void accumulateLightErrorSum(double* light_error_sum, double localreaction[2],
                                           double* oldreaction) {
  double x0 = oldreaction[0] - localreaction[0];
  double x1 = oldreaction[1] - localreaction[1];
  *light_error_sum += x0 * x0 + x1 * x1;
}
static double calculateLightError(double light_error_sum, unsigned int nc, double* reaction,
                                  double* norm_r) {
  double error = sqrt(light_error_sum);
  *norm_r = cblas_dnrm2((int32_t)nc * 2, reaction, 1);
  if (fabs(*norm_r) > DBL_EPSILON) error /= (*norm_r);
  return error;
}
static int determine_convergence(double error, double tolerance, unsigned int iter,
                                 SolverOptions* options) {
  int has_not_converged = 1;
  if (error < tolerance) {
    has_not_converged = 0;
    numerics_printf(
        "-- FC2D - NSGS - Iteration %i "
        "Residual = %14.7e < %7.3e\n",
        iter, error, tolerance);
  } else {
    numerics_printf(
        "-- FC2D - NSGS - Iteration %i "
        "Residual = %14.7e > %7.3e\n",
        iter, error, tolerance);
  }
  return has_not_converged;
}

static double calculateFullErrorFinal(FrictionContactProblem* problem, SolverOptions* options,
                                      double* reaction, double* velocity, double tolerance,
                                      double norm_q) {
  double absolute_error;
  /* (*computeError)(problem, reaction , velocity, tolerance, */
  /*                 options, norm_q, &absolute_error); */

  fc2d_compute_error(problem, reaction, velocity, tolerance, norm_q, &absolute_error);

  if (verbose > 0) {
    if (absolute_error > options->dparam[SICONOS_DPARAM_TOL]) {
      numerics_printf(
          "-- FC2D - NSGS - Warning absolute "
          "Residual = %14.7e is larger than required precision = %14.7e",
          absolute_error, options->dparam[SICONOS_DPARAM_TOL]);
    } else {
      numerics_printf(
          "-- FC2D - NSGS - absolute "
          "Residual = %14.7e is smaller than required precision = %14.7e",
          absolute_error, options->dparam[SICONOS_DPARAM_TOL]);
    }
  }
  return absolute_error;
}

static int determine_convergence_with_full_final(FrictionContactProblem* problem,
                                                 SolverOptions* options, double* reaction,
                                                 double* velocity, double* tolerance,
                                                 double norm_q, double error,
                                                 unsigned int iter) {
  int has_not_converged = 1;
  if (error < *tolerance) {
    has_not_converged = 0;
    numerics_printf(
        "-- FC2D - NSGS - Iteration %i "
        "Residual = %14.7e < %7.3e",
        iter, error, *tolerance);

    double absolute_error = calculateFullErrorFinal(
        problem, options, reaction, velocity, options->dparam[SICONOS_DPARAM_TOL], norm_q);
    if (absolute_error > options->dparam[SICONOS_DPARAM_TOL]) {
      *tolerance = error / absolute_error * options->dparam[SICONOS_DPARAM_TOL];
      assert(*tolerance > 0.0 && "tolerance has to be positive");
      /* if (*tolerance < DBL_EPSILON) */
      /* { */
      /*   numerics_warning("determine_convergence_with_full_fina", "We try to set a very smal
       * tolerance"); */
      /*   *tolerance = DBL_EPSILON; */
      /* } */
      numerics_printf(
          "-- FC2D - NSGS - We modify the required incremental precision to reach accuracy to "
          "%e",
          *tolerance);
      has_not_converged = 1;
    } else {
      numerics_printf(
          "-- FC2D - NSGS - The incremental precision is sufficient to reach accuracy to %e",
          *tolerance);
    }

  } else {
    numerics_printf(
        "-- FC2D - NSGS - Iteration %i "
        "Residual = %14.7e > %7.3e",
        iter, error, *tolerance);
  }
  return has_not_converged;
}
static double* fc2d_nsgs_compute_local_problem_determinant(
    SparseBlockCoordinateMatrix* diagonal_blocks) {
  double* diagonal_block_determinant =
      (double*)calloc(sizeof(double), diagonal_blocks->blocknumber0);
  for (unsigned int contact = 0; contact < diagonal_blocks->blocknumber0; ++contact) {
    double* block = diagonal_blocks->block[contact];
    diagonal_block_determinant[contact] = block[0] * block[3] - block[1] * block[2];
    if (diagonal_block_determinant[contact] < DBL_EPSILON) {
      free(diagonal_block_determinant);
      return NULL;
    }
  }
  return diagonal_block_determinant;
}

static inline void fc2d_nsgs_local_solve(double* W, double D, double* q, double mu,
                                         double* P) {
  /* | Wnn Wnt |
     | Wtn Wtt | */

#define Wnn W[0]
#define Wtn W[1]
#define Wnt W[2]
#define Wtt W[3]

  if (q[0] > 0) {
    P[0] = 0;
    P[1] = 0;
  } else {
    /* solve WP + q = 0  */

    P[0] = -(Wtt * q[0] - Wnt * q[1]) / D;
    P[1] = -(-Wtn * q[0] + Wnn * q[1]) / D;

    double muPn = mu * P[0];

    if (fabs(P[1]) > muPn)
    /* outside cone */
    {
      if (P[1] + muPn < 0) {
        P[0] = -q[0] / (Wnn - mu * Wnt);
        P[1] = -mu * P[0];
      } else {
        P[0] = -q[0] / (Wnn + mu * Wnt);
        P[1] = mu * P[0];
      }
    }
  }
#undef Wnn
#undef Wnt
#undef Wtn
#undef Wtt
}
static unsigned int* f2d_nsgs_allocate_freezing_contacts(FrictionContactProblem* problem,
                                                         SolverOptions* options) {
  unsigned int* fcontacts = 0;
  unsigned int nc = (unsigned int)problem->numberOfContacts;
  if (options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
    fcontacts = (unsigned int*)malloc(nc * sizeof(unsigned int));
    for (unsigned int i = 0; i < nc; ++i) {
      fcontacts[i] = 0;
    }
  }
  return fcontacts;
}

void fc2d_nsgs_permut(FrictionContactProblem* problem, double* z, double* w, int* info,
                      SolverOptions* options) {
  /* Notes:
     - we suppose that the trivial solution case has been checked before,
     and that all inputs differs from NULL since this function is
     supposed to be called from lcp_driver_global().
  */
  /* verbose=1; */
  /* Global Solver parameters*/
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = dparam[SICONOS_DPARAM_TOL];

  unsigned int nc = (unsigned int)problem->numberOfContacts;
  double norm_q = cblas_dnrm2((int32_t)nc * 2, problem->q, 1);
  double norm_r[] = {INFINITY};

  // Get coloring and permutation
  size_t n_colors = 0;
  size_t* sum_sizes = NULL;
  size_t* inv_permutation = (size_t*)malloc((size_t)nc * sizeof(size_t));
  color_graph_block_permut((size_t)problem->numberOfContacts, problem->M, &n_colors,
                           &sum_sizes, inv_permutation);

  free(sum_sizes);

  /* Local problem initialization */
  LinearComplementarityProblem* local_problem =
      (LinearComplementarityProblem*)malloc(sizeof(*local_problem));

  SparseBlockCoordinateMatrix* diagonal_blocks = fc3d_extract_diagonal_blocks(problem);

  double* diagonal_block_determinant =
      fc2d_nsgs_compute_local_problem_determinant(diagonal_blocks);
  /* verbose if problem */
  if (!diagonal_block_determinant) {
    /* Number of GS iterations */
    iparam[SICONOS_IPARAM_ITER_DONE] = 0;
    dparam[SICONOS_DPARAM_RESIDU] = INFINITY;
    numerics_printf("-- FC2D - NSGS - error: determinant diagonal block in W is zero \n");
    *info = 1;
    diagonal_blocks = fc3d_free_diagonal_blocks(problem, diagonal_blocks);
    free(diagonal_block_determinant);
    free(local_problem->q);
    free(local_problem->M);
    free(local_problem);
  }

  local_problem->M = NM_new();
  local_problem->M->storageType = NM_DENSE;
  local_problem->M->size0 = 2;
  local_problem->M->size1 = 2;

  local_problem->q = (double*)malloc(2 * sizeof(double));

  double localreaction[2];

  /*****  Gauss-Seidel iterations *****/
  int iter = 0;            /* Current iteration number */
  double error = INFINITY; /* Current error */
  int has_not_converged = 1;

  unsigned int* freeze_contacts = NULL;
  if (iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
    freeze_contacts = f2d_nsgs_allocate_freezing_contacts(problem, options);

    while ((iter < itermax) && has_not_converged) {
      ++iter;

      double light_error_sum = 0.0;
      double light_error_2 = 0.0;
      /* Loop over the rows of blocks in blmat */
      /* contact: current row (of blocks) number */
      unsigned int number_of_freezed_contact = 0;
      double tmp_criteria1 = tolerance * tolerance / (nc * nc * 10);
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
      for (unsigned int pos = 0, contact = 0; contact < nc; ++contact, ++pos, ++pos) {
        if (freeze_contacts[contact] > 0) {
          /* we skip freeze contacts */
          freeze_contacts[contact] -= 1;
          continue;
        }

        /* store  old reaction */
        localreaction[0] = z[pos];
        localreaction[1] = z[pos + 1];

        /* Local problem formalization */
        fc2d_nsgs_buildLocalProblem((int)contact, problem, diagonal_blocks, local_problem, z);

        /* Solve local problem */
        fc2d_nsgs_local_solve(local_problem->M->matrix0, diagonal_block_determinant[contact],
                              local_problem->q, problem->mu[contact], localreaction);

        light_error_2 = light_error_squared(localreaction, &z[pos]);
        light_error_sum += light_error_2;
        int relative_convergence_criteria =
            light_error_2 <= tmp_criteria1 * squared_norm(localreaction);
        int small_reaction_criteria = squared_norm(localreaction) <= tmp_criteria2;
        if ((relative_convergence_criteria || small_reaction_criteria) && iter >= 10) {
          /* we  freeze the contact for n iterations*/
          freeze_contacts[contact] =
              (unsigned int)iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT];
          DEBUG_EXPR(printf("first criteria : light_error_2*squared_norm(localreaction) <= "
                            "tolerance*tolerance/(nc*nc*10) ==> %e <= %e, bool =%i\n",
                            light_error_2 * squared_norm(localreaction),
                            tolerance * tolerance / (nc * nc * 10),
                            relative_convergence_criteria);
                     printf("second criteria :  squared_norm(localreaction) <=  (*norm_r* "
                            "*norm_r/(nc*nc))/1000. ==> %e <= %e, bool =%i \n",
                            squared_norm(localreaction), *norm_r * *norm_r / (nc * nc * 1000),
                            small_reaction_criteria);
                     printf("Contact % i is freezed for %i steps\n", contact,
                            iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT]););
        }
        /* reaction update */
        z[pos] = localreaction[0];
        z[pos + 1] = localreaction[1];

      }  // end for loop

      DEBUG_EXPR(int frozen_contact = 0;
                 for (unsigned int ii = 0; ii < nc; ++ii) if (freeze_contacts[ii] > 0)
                     frozen_contact++;
                 numerics_printf_verbose(1, "number of frozen contacts %i at iter : %i",
                                         frozen_contact, iter););

      /* error evaluation */
      if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
          SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT) {
        error = calculateLightError(light_error_sum, nc, z, norm_r);
        has_not_converged =
            determine_convergence(error, tolerance, (unsigned int)iter, options);
      } else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                 SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
        error = calculateLightError(light_error_sum, nc, z, norm_r);
        has_not_converged = determine_convergence_with_full_final(
            problem, options, z, w, &tolerance, norm_q, error, (unsigned int)iter);
      }
    }  // end while loop

  } else {
    while ((iter < itermax) && has_not_converged) {
      double light_error_sum = 0.0;
      unsigned int pos;
      unsigned int contact;
      /* Loop over the rows of blocks in blmat */
      for (unsigned int i = 0; i < nc; i++) {
        contact = (unsigned int)inv_permutation[i];
        pos = 2 * contact;

        /* store  old reaction */
        localreaction[0] = z[pos];
        localreaction[1] = z[pos + 1];

        /* Local problem formalization */
        fc2d_nsgs_buildLocalProblem((int)contact, problem, diagonal_blocks, local_problem, z);

        /* Solve local problem */
        fc2d_nsgs_local_solve(local_problem->M->matrix0, diagonal_block_determinant[contact],
                              local_problem->q, problem->mu[contact], localreaction);

        if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT ||
            iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
          accumulateLightErrorSum(&light_error_sum, localreaction, &z[pos]);

        z[pos] = localreaction[0];
        z[pos + 1] = localreaction[1];

      }  // end for loop

      /*  error evaluation */
      if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
          SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT) {
        error = calculateLightError(light_error_sum, nc, z, norm_r);
        has_not_converged =
            determine_convergence(error, tolerance, (unsigned int)iter, options);
      } else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                 SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
        error = calculateLightError(light_error_sum, nc, z, norm_r);
        has_not_converged = determine_convergence_with_full_final(
            problem, options, z, w, &tolerance, norm_q, error, (unsigned int)iter);
      }
      ++iter;
    }  // end while loop
  }

  /* Full criterium */
  if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
      SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
    error = calculateFullErrorFinal(problem, options, z, w, tolerance, norm_q);

    has_not_converged =
        determine_convergence(error, dparam[SICONOS_DPARAM_TOL], (unsigned int)iter, options);
  }

  // numerics_printf("Siconos Numerics : problem size=%d, nb iterations=%d, error=%g\n",
  //          blmat->blocknumber0,
  //          iter,
  //          error);

  *info = has_not_converged;

  /* Number of GS iterations */
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

  /* Resulting error */
  dparam[SICONOS_DPARAM_RESIDU] = error;

  if (freeze_contacts) free(freeze_contacts);
  diagonal_blocks = fc3d_free_diagonal_blocks(problem, diagonal_blocks);
  free(diagonal_block_determinant);
  free(local_problem->q);
  free(local_problem->M);
  free(local_problem);
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers FC2D_NSGS in the global solver registry, enabling:
 * - Dynamic solver lookup by ID
 * - Runtime solver introspection
 * - Elimination of giant switch statements in drivers
 */

static int fc2d_nsgs_init_wrap(void* problem, SolverOptions* options) {
  fc2d_nsgs_set_default(options);
  return NUMERICS_OK;
}

static int fc2d_nsgs_solve_wrap(void* problem, double* reaction, double* velocity,
                                SolverOptions* options) {
  int info = NUMERICS_OK;
  fc2d_nsgs((FrictionContactProblem*)problem, reaction, velocity, &info, options);
  return info;
}

static void fc2d_nsgs_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_FRICTION_2D_NSGS_PERMUT, "SICONOS_FRICTION_2D_NSGS_PERMUT",
                "Test solver: sequential FC2D_NSGS with the parallel contact order",
                fc2d_nsgs_init_wrap, fc2d_nsgs_solve_wrap, fc2d_nsgs_free_wrap,
                NULL,                  /* error function */
                fc2d_nsgs_set_default, /* set_default */
                1000,                  /* default_max_iter */
                1e-4,                  /* default_tol */
                0)                     /* is_local_solver */