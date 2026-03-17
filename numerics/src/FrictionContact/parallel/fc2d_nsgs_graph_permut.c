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
#include <string.h>  // for memcpy

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
#include "op3x3.h"

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

/* Only support SBM format for now */
static double* fc2d_extract_diagonal_blocks(FrictionContactProblem* problem) {
  unsigned int nc = (unsigned int)problem->numberOfContacts;
  double* sbcm = (double*)calloc(4 * nc, sizeof(double));
  double* block;
  int diagPos;

  for (unsigned int i = 0; i < nc; ++i) {
    diagPos = (int)SBM_diagonal_block_index(problem->M->matrix1, i);
    block = problem->M->matrix1->block[diagPos];
    for (unsigned int j = 0; j < 4; j++) sbcm[i * 4 + j] = block[j];
  }

  return sbcm;
}
static void* fc2d_free_diagonal_blocks(double* sbcm) {
  free(sbcm);
  return NULL;
}

static void fc2d_nsgs_buildLocalProblem_parallel(
    unsigned int contact, FrictionContactProblem* problem, double* blocks_contiguous,
    double* diagonal_blocks, size_t* index1_data, size_t* index2_data,
    LinearComplementarityProblem* local_problem, double* reaction) {
  // NM_extract_diag_block2(problem->M, contact, &local_problem->M->matrix0);
  local_problem->M->matrix0 = &diagonal_blocks[4 * contact];

  local_problem->M->size0 = 2;  // Necessary ?
  local_problem->M->size1 = 2;

  local_problem->q[0] = problem->q[contact * 2];
  local_problem->q[1] = problem->q[contact * 2 + 1];

  size_t colNumber;

  for (size_t blockNum = index1_data[contact]; blockNum < index1_data[contact + 1];
       ++blockNum) {
    colNumber = index2_data[blockNum];
    if (colNumber != contact) {
      mvp2x2(&blocks_contiguous[blockNum * 4], &reaction[2 * colNumber], local_problem->q);
    }
  }

  /* NM_row_prod_no_diag2_parallel(2 * problem->numberOfContacts, contact, 2 * contact,
     problem->M, reaction, local_problem->q, false); */
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

static double* fc2d_nsgs_compute_local_problem_determinant(unsigned int nc,
                                                           double* diagonal_blocks) {
  double* diagonal_block_determinant = (double*)calloc(sizeof(double), nc);
  double* block;

  for (unsigned int i = 0; i < nc; ++i) {
    block = &diagonal_blocks[4 * i];
    diagonal_block_determinant[i] = block[0] * block[3] - block[1] * block[2];
    if (diagonal_block_determinant[i] < DBL_EPSILON) {
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

void fc2d_nsgs_graph_permut(FrictionContactProblem* problem, double* z, double* w, int* info,
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
  double norm_r = 0.0;

  /* COLORING AND PERMUTATION */

  /* Coloring */
  size_t n_colors = 0;
  size_t* sum_sizes = NULL;
  size_t* inv_permutation = (size_t*)malloc((size_t)nc * sizeof(size_t));
  color_graph_block_permut((size_t)problem->numberOfContacts, problem->M, &n_colors,
                           &sum_sizes, inv_permutation);

  SparseBlockStructuredMatrix* SBM_problem = NULL;

  /* Create SBM if it does not exist */
  if (!problem->M->matrix1) {
    SBM_problem = SBM_new();
    switch (problem->M->storageType) {
      // Convert DENSE -> SBM
      case NM_DENSE: {
        SBM_from_dense(problem->dimension, (size_t)problem->M->size1,
                       (size_t)problem->M->size0, problem->M->matrix0, SBM_problem);
        /* printf("Number of blocks: %d\n", SBM_problem->nbblocks);
        for (int b = 0; b < SBM_problem->nbblocks; b++) {
          printf("%e %e %e %e\n", SBM_problem->block[b][0], SBM_problem->block[b][1],
                 SBM_problem->block[b][2], SBM_problem->block[b][3]);
        } */
        break;
      }
      // Convert SPARSE -> SBM
      case NM_SPARSE: {
        // Get Csparse matrix
        CSparseMatrix* sparse = NULL;
        if (problem->M->matrix2->origin == NSM_CSR) {
          sparse = NM_csr(problem->M);
        } else {
          sparse = NM_csc_trans(problem->M);
        }

        SBM_from_csparse_2(problem->dimension, sparse, SBM_problem);

        break;
      }
      // Add this cases to avoid warning
      case NM_UNKNOWN:
      case NM_SPARSE_BLOCK: {
        break;
      }
    }
    problem->M->matrix1 = SBM_problem;
  }

  /* Switch storageType to SBM temporarily */
  unsigned int old_storageType = problem->M->storageType;
  problem->M->storageType = NM_SPARSE_BLOCK;

  /* Permutate rows and columns of SBM */
  /* Can  do better? In place stuff? */
  SparseBlockStructuredMatrix* SBM_col_permuted = SBM_new();
  SparseBlockStructuredMatrix* SBM_permuted = SBM_new();
  unsigned int* rowIndex = (unsigned int*)malloc(nc * sizeof(unsigned int));
  for (unsigned int i = 0; i < nc; i++) rowIndex[inv_permutation[i]] = i;

  SBM_column_permutation(rowIndex, problem->M->matrix1, SBM_col_permuted);
  SBM_row_permutation_copy(inv_permutation, SBM_col_permuted, SBM_permuted);

  free(rowIndex);

  SparseBlockStructuredMatrix* old_matrix1 = problem->M->matrix1;
  problem->M->matrix1 = SBM_permuted;

  /* Store all blocks in contiguous array */
  unsigned int nbblocks = problem->M->matrix1->nbblocks;
  double* blocks_contiguous = (double*)malloc(nbblocks * 4 * sizeof(double));
  double* current_block;
  for (unsigned int blockNum = 0; blockNum < nbblocks; blockNum++) {
    current_block = problem->M->matrix1->block[blockNum];
    for (unsigned int j = 0; j < 4; j++) {
      blocks_contiguous[blockNum * 4 + j] = current_block[j];
    }
  }

  double* mu_permuted = (double*)malloc(nc * sizeof(double));
  double* q_permuted = (double*)malloc(nc * 2 * sizeof(double));
  for (unsigned int i = 0; i < nc; i++) {
    q_permuted[2 * i] = problem->q[2 * inv_permutation[i]];
    q_permuted[2 * i + 1] = problem->q[2 * inv_permutation[i] + 1];
    mu_permuted[i] = problem->mu[inv_permutation[i]];
  }

  double* old_q = problem->q;
  double* old_mu = problem->mu;
  problem->q = q_permuted;
  problem->mu = mu_permuted;

  /* Get diagonal blocks and determinants */
  double* diagonal_blocks = fc2d_extract_diagonal_blocks(problem);
  double* diagonal_block_determinant =
      fc2d_nsgs_compute_local_problem_determinant(nc, diagonal_blocks);

  /* Local problem initialization */
  LinearComplementarityProblem* local_problem =
      (LinearComplementarityProblem*)malloc(sizeof(LinearComplementarityProblem));
  local_problem->M = NM_new();
  local_problem->M->storageType = NM_DENSE;
  local_problem->M->size0 = 2;
  local_problem->M->size1 = 2;
  local_problem->q = (double*)malloc(2 * sizeof(double));

  /*****  Gauss-Seidel iterations *****/
  int iter = 0;            /* Current iteration number */
  double error = INFINITY; /* Current error */
  int has_not_converged = 1;

  size_t* index1_data = NULL;
  size_t* index2_data = NULL;

  unsigned int* freeze_contacts = NULL;
  // FREEZING CONTACTS
  if (iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
    unsigned int pos;
    double light_error_sum;
    double light_error_2;
    double localreaction[2];
    unsigned int number_of_freezed_contact = 0;
    double tmp_criteria1, tmp_criteria2;

    freeze_contacts = f2d_nsgs_allocate_freezing_contacts(problem, options);

#pragma omp parallel default(none)                                                      \
    private(pos, local_problem, localreaction, light_error_2, index1_data, index2_data) \
    shared(problem, diagonal_blocks, diagonal_block_determinant, z, sum_sizes, iter)    \
    shared(iparam, light_error_sum, n_colors, norm_r, nc, error, options, tolerance,    \
               has_not_converged, norm_q, w, itermax)                                   \
    shared(tmp_criteria1, tmp_criteria2, freeze_contacts, number_of_freezed_contact,    \
               blocks_contiguous)
    {
      local_problem =
          (LinearComplementarityProblem*)malloc(sizeof(LinearComplementarityProblem));
      local_problem->M = NM_new();
      local_problem->M->storageType = NM_DENSE;
      local_problem->M->size0 = 2;
      local_problem->M->size1 = 2;
      local_problem->q = (double*)malloc(2 * sizeof(double));

      index1_data = (size_t*)malloc((nc + 1) * sizeof(size_t));
      index2_data = (size_t*)malloc(problem->M->matrix1->nbblocks * sizeof(size_t));

      memcpy(index1_data, problem->M->matrix1->index1_data, (nc + 1) * sizeof(size_t));
      memcpy(index2_data, problem->M->matrix1->index2_data,
             problem->M->matrix1->nbblocks * sizeof(size_t));

      while ((iter < itermax) && has_not_converged) {
        // light_error_sum = 0.0;
        light_error_2 = 0.0;
        /* Loop over the rows of blocks in blmat */
        /* contact: current row (of blocks) number */

#pragma omp single
        {
          number_of_freezed_contact = 0;
          tmp_criteria1 = tolerance * tolerance / (nc * nc * 10);
          tmp_criteria2 = norm_r * norm_r / (nc * nc * 1000);
          if (iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
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
#pragma omp for reduction(+ : light_error_sum)
          for (unsigned int permuted_contact = (unsigned int)sum_sizes[color];
               permuted_contact < sum_sizes[color + 1]; permuted_contact++) {
            if (freeze_contacts[permuted_contact] > 0) {
              /* we skip freeze contacts */
              freeze_contacts[permuted_contact] -= 1;
              continue;
            }

            pos = 2 * permuted_contact;
            localreaction[0] = z[pos];
            localreaction[1] = z[pos + 1];

            /* Local problem formalization */
            fc2d_nsgs_buildLocalProblem_parallel(permuted_contact, problem, blocks_contiguous,
                                                 diagonal_blocks, index1_data, index2_data,
                                                 local_problem, z);

            /* Solve local problem */
            fc2d_nsgs_local_solve(
                local_problem->M->matrix0, diagonal_block_determinant[permuted_contact],
                local_problem->q, problem->mu[permuted_contact], localreaction);

            light_error_2 = light_error_squared(localreaction, &z[pos]);

            // #pragma omp atomic update
            light_error_sum += light_error_2;

            int relative_convergence_criteria =
                light_error_2 <= tmp_criteria1 * squared_norm(localreaction);
            int small_reaction_criteria = squared_norm(localreaction) <= tmp_criteria2;
            if ((relative_convergence_criteria || small_reaction_criteria) && iter >= 10) {
              /* we  freeze the contact for n iterations*/
              freeze_contacts[permuted_contact] =
                  (unsigned int)iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT];
              DEBUG_EXPR(
                  printf("first criteria : light_error_2*squared_norm(localreaction) <= "
                         "tolerance*tolerance/(nc*nc*10) ==> %e <= %e, bool =%i\n",
                         light_error_2 * squared_norm(localreaction),
                         tolerance * tolerance / (nc * nc * 10),
                         relative_convergence_criteria);
                  printf("second criteria :  squared_norm(localreaction) <=  (*norm_r* "
                         "*norm_r/(nc*nc))/1000. ==> %e <= %e, bool =%i \n",
                         squared_norm(localreaction), norm_r * norm_r / (nc * nc * 1000),
                         small_reaction_criteria);
                  printf("Contact % i is freezed for %i steps\n", permuted_contact,
                         iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT]););
            }
            /* reaction update */
            z[pos] = localreaction[0];
            z[pos + 1] = localreaction[1];
          }
        }  // end for loop

#pragma omp single
        {
          DEBUG_EXPR(int frozen_contact = 0;
                     for (unsigned int ii = 0; ii < nc; ++ii) if (freeze_contacts[ii] > 0)
                         frozen_contact++;
                     numerics_printf_verbose(1, "number of frozen contacts %i at iter : %i",
                                             frozen_contact, iter););

          /* error evaluation */
          if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
              SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT) {
            error = calculateLightError(light_error_sum, nc, z, &norm_r);
            has_not_converged =
                determine_convergence(error, tolerance, (unsigned int)iter, options);
          } else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                     SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
            error = calculateLightError(light_error_sum, nc, z, &norm_r);
            has_not_converged = determine_convergence_with_full_final(
                problem, options, z, w, &tolerance, norm_q, error, (unsigned int)iter);
            problem->M->storageType = NM_SPARSE_BLOCK;
          }
          light_error_sum = 0.;
          ++iter;
        }
      }  // end while loop
      free(local_problem->q);
      free(local_problem->M);
      free(local_problem);
    }  // end parallel region

    /***********************/
    /* NO FREEZING CONTACT */
    /***********************/
  } else {
    unsigned int pos;
    double light_error_sum = 0.;
    double localreaction[2];

#pragma omp parallel default(none)                                                   \
    private(pos, local_problem, localreaction, index1_data, index2_data)             \
    shared(problem, diagonal_blocks, diagonal_block_determinant, z, sum_sizes, iter, \
               blocks_contiguous)                                                    \
    shared(iparam, light_error_sum, n_colors, norm_r, nc, error, options, tolerance, \
               has_not_converged, norm_q, w, itermax)
    {
      local_problem = (LinearComplementarityProblem*)malloc(sizeof(*local_problem));
      local_problem->M = NM_new();
      local_problem->M->storageType = NM_DENSE;
      local_problem->M->size0 = 2;
      local_problem->M->size1 = 2;
      local_problem->q = (double*)malloc(2 * sizeof(double));

      index1_data = (size_t*)malloc((nc + 1) * sizeof(size_t));
      index2_data = (size_t*)malloc(problem->M->matrix1->nbblocks * sizeof(size_t));

      memcpy(index1_data, problem->M->matrix1->index1_data, (nc + 1) * sizeof(size_t));
      memcpy(index2_data, problem->M->matrix1->index2_data,
             problem->M->matrix1->nbblocks * sizeof(size_t));

      double* blocks_contiguous_local =
          (double*)malloc(problem->M->matrix1->nbblocks * 4 * sizeof(double));
      memcpy(blocks_contiguous_local, blocks_contiguous,
             problem->M->matrix1->nbblocks * 4 * sizeof(double));
      double* diagonal_blocks_local = fc2d_extract_diagonal_blocks(problem);

      double local_light_error_sum = 0.0;
      double local_norm_r = 0.0;

      // unsigned int n_thread = omp_get_thread_num();

      while ((iter < itermax) && has_not_converged) {
        for (size_t color = 0; color < n_colors; color++) {
#pragma omp for
          for (unsigned int permuted_contact = (unsigned int)sum_sizes[color];
               permuted_contact < sum_sizes[color + 1]; permuted_contact++) {
            pos = 2 * permuted_contact;
            localreaction[0] = z[pos];
            localreaction[1] = z[pos + 1];

            /* Local problem formalization */
            fc2d_nsgs_buildLocalProblem_parallel(
                permuted_contact, problem, blocks_contiguous_local, diagonal_blocks_local,
                index1_data, index2_data, local_problem, z);

            /* Solve local problem */
            fc2d_nsgs_local_solve(
                local_problem->M->matrix0, diagonal_block_determinant[permuted_contact],
                local_problem->q, problem->mu[permuted_contact], localreaction);

            local_light_error_sum += light_error_squared(localreaction, &z[pos]);
            local_norm_r += squared_norm(localreaction);

            /* Update z */
            z[pos] = localreaction[0];
            z[pos + 1] = localreaction[1];
          }

        }  // end for loop

        // manual_reduction[2 * n_thread] = local_light_error_sum;
        // manual_reduction[2 * n_thread + 1] = local_norm_r;

#pragma omp atomic
        light_error_sum += local_light_error_sum;
#pragma omp atomic
        norm_r += local_norm_r;

#pragma omp barrier

        local_light_error_sum = 0.0;
        local_norm_r = 0.0;

/*
#pragma omp atomic
light_error_sum += local_light_error_sum;
#pragma omp atomic
norm_r += local_norm_r;

local_light_error_sum = 0.0;
local_norm_r = 0.0;
*/

/*  error evaluation */
#pragma omp single
        {
          // error = calculateLightError(light_error_sum, nc, z, norm_r);
          error = sqrt(light_error_sum);
          norm_r = sqrt(norm_r);
          if (fabs(norm_r) > DBL_EPSILON) error /= norm_r;

          if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
              SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT) {
            has_not_converged =
                determine_convergence(error, tolerance, (unsigned int)iter, options);
          } else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                     SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
            has_not_converged = determine_convergence_with_full_final(
                problem, options, z, w, &tolerance, norm_q, error, (unsigned int)iter);
            problem->M->storageType = NM_SPARSE_BLOCK;
          }

          light_error_sum = 0.;
          norm_r = 0.0;
          ++iter;
        }

      }  // end while loop
      free(local_problem->q);
      free(local_problem->M);
      free(local_problem);
      free(index1_data);
      free(index2_data);
      free(blocks_contiguous_local);
      free(diagonal_blocks_local);
    }
  }

  /* Full criterium */
  if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
      SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
    error = calculateFullErrorFinal(problem, options, z, w, tolerance, norm_q);
    problem->M->storageType = NM_SPARSE_BLOCK;

    has_not_converged =
        determine_convergence(error, dparam[SICONOS_DPARAM_TOL], (unsigned int)iter, options);
  }

  // numerics_printf("Siconos Numerics : problem size=%d, nb iterations=%d, error=%g\n",
  //          blmat->blocknumber0,
  //          iter,
  //          error);

  *info = has_not_converged;

  /* Permutate z and w back */
  double* z_permut = (double*)malloc(2 * nc * sizeof(double));
  double* w_permut = (double*)malloc(2 * nc * sizeof(double));
  for (unsigned int i = 0; i < nc; i++) {
    z_permut[2 * inv_permutation[i]] = z[2 * i];
    z_permut[2 * inv_permutation[i] + 1] = z[2 * i + 1];
    w_permut[2 * inv_permutation[i]] = w[2 * i];
    w_permut[2 * inv_permutation[i] + 1] = w[2 * i + 1];
  }
  memcpy(z, z_permut, 2 * nc * sizeof(double));
  memcpy(w, w_permut, 2 * nc * sizeof(double));
  free(z_permut);
  free(w_permut);

  /* Number of GS iterations */
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

  /* Resulting error */
  dparam[SICONOS_DPARAM_RESIDU] = error;

  if (freeze_contacts) free(freeze_contacts);
  diagonal_blocks = fc2d_free_diagonal_blocks(diagonal_blocks);
  free(diagonal_block_determinant);
  free(local_problem->q);
  free(local_problem->M);
  free(local_problem);

  free(sum_sizes);
  free(inv_permutation);

  /* Restore problem before permutation */
  problem->M->matrix1 = old_matrix1;
  problem->q = old_q;
  problem->mu = old_mu;
  problem->M->storageType = old_storageType;

  free(blocks_contiguous);

  if (SBM_problem != NULL) {
    SBM_clear(SBM_problem);
    free(SBM_problem);
    SBM_problem = NULL;
    problem->M->matrix1 = NULL;
  }
  // free(SBM_problem);

  SBMfree(SBM_col_permuted, 0);  // do not free blocks on this one
  SBM_clear(SBM_permuted);       // free blocks because they were copied
  free(SBM_col_permuted);
  free(SBM_permuted);
  free(q_permuted);
  free(mu_permuted);

  /* DO NOT FORGET TO FREE THE REST */
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

REGISTER_SOLVER(SICONOS_FRICTION_2D_NSGS_GRAPH_PERMUT, "SICONOS_FRICTION_2D_NSGS_GRAPH_PERMUT",
                "Parallel version of FC2D_NSGS (OpenMP, with permutation)",
                fc2d_nsgs_init_wrap, fc2d_nsgs_solve_wrap, fc2d_nsgs_free_wrap,
                NULL,                  /* error function */
                fc2d_nsgs_set_default, /* set_default */
                1000,                  /* default_max_iter */
                1e-4,                  /* default_tol */
                0)                     /* is_local_solver */