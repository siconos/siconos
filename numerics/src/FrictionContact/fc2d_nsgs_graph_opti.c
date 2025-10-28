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
#include "Friction_cst.h"                  // for SICONOS_FRICTION_3D_IPARAM...
#include "LCP_Solvers.h"                   // for lcp_nsgs_SBM_buildLocalPro...
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NumericsFwd.h"                   // for SolverOptions, LinearCompl...
#include "NumericsMatrix.h"                // for NumericsMatrix, RawNumeric...
#include "SiconosBlas.h"                   // for cblas_dnrm2
#include "SolverOptions.h"                 // for SolverOptions, SICONOS_DPA...
#include "SparseBlockMatrix.h"             // for SparseBlockStructuredMatrix
#include "fc2d_Solvers.h"                  // for fc2d_nsgs_sbm, fc2d_spa...
#include "fc2d_compute_error.h"            // for fc2d_compute_error
#include "numerics_verbose.h"              // for numerics_printf, verbose
#include "graph_tools.h"

#include "op3x3.h"

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES 1 */
#include "siconos_debug.h"  // for DEBUG_BEGIN, DEBUG_END
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

#define SGN(x) ((x) < 0 ? -1 : (x) > 0 ? 1 : 0)

static SparseBlockCoordinateMatrix *fc3d_extract_diagonal_blocks(
    FrictionContactProblem *problem) {
  unsigned int nc = problem->numberOfContacts;

  SparseBlockCoordinateMatrix *sbcm = SBCM_new();
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
  sbcm->block = (double **)malloc(sizeof(double *) * nbblocks);
  //  sbcm->blocksize0 = (unsigned int *) malloc(sizeof(unsigned int) * nbblocks);
  //  sbcm->blocksize1 = (unsigned int *) malloc(sizeof(unsigned int) * nbblocks);

  for (unsigned int contact = 0; contact < nc; ++contact) {
    //    sbcm->blocksize0[contact]=2;
    //    sbcm->blocksize1[contact]=2;
    if (problem->M->storageType != NM_SPARSE_BLOCK)
      sbcm->block[contact] = (double *)calloc(sizeof(double *), 4);
    NM_extract_diag_block2(problem->M, contact, &sbcm->block[contact]);
  }

  return sbcm;
}
static SparseBlockCoordinateMatrix *fc3d_free_diagonal_blocks(
    FrictionContactProblem *problem, SparseBlockCoordinateMatrix *sbcm) {
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

static void fc2d_nsgs_buildLocalProblem(int contact, FrictionContactProblem *problem,
                                        SparseBlockCoordinateMatrix *diagonal_blocks,
                                        LinearComplementarityProblem *local_problem,
                                        double *reaction) {
  // NM_extract_diag_block2(problem->M, contact, &local_problem->M->matrix0);
  local_problem->M->matrix0 = diagonal_blocks->block[contact];

  local_problem->M->size0 = 2;  // Necessary ?
  local_problem->M->size1 = 2;

  local_problem->q[0] = problem->q[contact * 2];
  local_problem->q[1] = problem->q[contact * 2 + 1];
  NM_row_prod_no_diag2(2 * problem->numberOfContacts, contact, 2 * contact, problem->M,
                       reaction, local_problem->q, false);

  DEBUG_EXPR(NM_display(local_problem->M););
  DEBUG_EXPR(NV_display(local_problem->q, 2););
}

static void fc2d_nsgs_buildLocalProblem_parallel(int contact, FrictionContactProblem *problem,
                                                 SparseBlockCoordinateMatrix *diagonal_blocks,
                                                 LinearComplementarityProblem *local_problem,
                                                 double *reaction) {
  // NM_extract_diag_block2(problem->M, contact, &local_problem->M->matrix0);
  local_problem->M->matrix0 = diagonal_blocks->block[contact];

  local_problem->M->size0 = 2;  // Necessary ?
  local_problem->M->size1 = 2;

  local_problem->q[0] = problem->q[contact * 2];
  local_problem->q[1] = problem->q[contact * 2 + 1];
  NM_row_prod_no_diag2_parallel(2 * problem->numberOfContacts, contact, 2 * contact, problem->M,
                                reaction, local_problem->q, false);

  DEBUG_EXPR(NM_display(local_problem->M););
  DEBUG_EXPR(NV_display(local_problem->q, 2););
}

static void shuffle(unsigned int size, unsigned int *randnum)  // size is the given range
{
  unsigned int swap, randindex;
  for (unsigned i = 0; i < size; ++i) {
    swap = randnum[i];
    randindex = rand() % size;
    randnum[i] = randnum[randindex];
    randnum[randindex] = swap;
  }
}

static inline double light_error_squared(double localreaction[2], double *oldreaction) {
  double x0 = oldreaction[0] - localreaction[0];
  double x1 = oldreaction[1] - localreaction[1];
  return x0 * x0 + x1 * x1;
}
static inline double squared_norm(double localreaction[2]) {
  return (localreaction[0] * localreaction[0] + localreaction[1] * localreaction[1]);
}

static inline void accumulateLightErrorSum(double *light_error_sum, double localreaction[2],
                                           double *oldreaction) {
  double x0 = oldreaction[0] - localreaction[0];
  double x1 = oldreaction[1] - localreaction[1];
  *light_error_sum += x0 * x0 + x1 * x1;
}
static double calculateLightError(double light_error_sum, unsigned int nc, double *reaction,
                                  double *norm_r) {
  double error = sqrt(light_error_sum);
  *norm_r = cblas_dnrm2(nc * 2, reaction, 1);
  if (fabs(*norm_r) > DBL_EPSILON) error /= (*norm_r);
  return error;
}
static int determine_convergence(double error, double tolerance, unsigned int iter,
                                 SolverOptions *options) {
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

static double calculateFullErrorFinal(FrictionContactProblem *problem, SolverOptions *options,
                                      double *reaction, double *velocity, double tolerance,
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

static int determine_convergence_with_full_final(FrictionContactProblem *problem,
                                                 SolverOptions *options, double *reaction,
                                                 double *velocity, double *tolerance,
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
static double *fc2d_nsgs_compute_local_problem_determinant(
    SparseBlockCoordinateMatrix *diagonal_blocks) {
  double *diagonal_block_determinant =
      (double *)calloc(sizeof(double), diagonal_blocks->blocknumber0);
  for (unsigned int contact = 0; contact < diagonal_blocks->blocknumber0; ++contact) {
    double *block = diagonal_blocks->block[contact];
    diagonal_block_determinant[contact] = block[0] * block[3] - block[1] * block[2];
    if (diagonal_block_determinant[contact] < DBL_EPSILON) {
      free(diagonal_block_determinant);
      return NULL;
    }
  }
  return diagonal_block_determinant;
}

static inline void fc2d_nsgs_local_solve(double *W, double D, double *q, double mu,
                                         double *P) {
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
static unsigned int *f2d_nsgs_allocate_freezing_contacts(FrictionContactProblem *problem,
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

void fc2d_nsgs_graph_opti(FrictionContactProblem *problem, double *z, double *w, int *info,
                          SolverOptions *options) {
  /* Notes:
     - we suppose that the trivial solution case has been checked before,
     and that all inputs differs from NULL since this function is
     supposed to be called from lcp_driver_global().
  */
  /* verbose=1; */
  /* Global Solver parameters*/
  int *iparam = options->iparam;
  double *dparam = options->dparam;

  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = dparam[SICONOS_DPARAM_TOL];

  unsigned int nc = problem->numberOfContacts;
  double norm_q = cblas_dnrm2(nc * 2, problem->q, 1);
  double norm_r[] = {INFINITY};

  /* Local problem initialization */
  LinearComplementarityProblem *local_problem =
      (LinearComplementarityProblem *)malloc(sizeof(*local_problem));

  SparseBlockCoordinateMatrix *diagonal_blocks = fc3d_extract_diagonal_blocks(problem);

  double *diagonal_block_determinant =
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

  local_problem->q = (double *)malloc(2 * sizeof(double));

  /* Coloring */
  size_t n_colors = 0;
  size_t *partition_size = NULL;
  size_t **partitions = NULL;
  color_graph_block(problem->numberOfContacts, problem->M, &n_colors, &partition_size, &partitions);

  /*****  Gauss-Seidel iterations *****/
  int iter = 0;            /* Current iteration number */
  double error = INFINITY; /* Current error */
  int has_not_converged = 1;

  unsigned int *freeze_contacts = NULL;
  // FREEZING CONTACTS
  if (iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
    unsigned int contact;
    unsigned int pos;
    double light_error_sum;
    double light_error_2;
    double localreaction[2];
    unsigned int number_of_freezed_contact = 0;
    double tmp_criteria1, tmp_criteria2;

    freeze_contacts = f2d_nsgs_allocate_freezing_contacts(problem, options);

    #pragma omp parallel default(none) \
                         private(contact, pos, local_problem, localreaction, light_error_2) \
                         shared(problem, diagonal_blocks, diagonal_block_determinant, z, partition_size, partitions, iter) \
                         shared(iparam, light_error_sum, n_colors, norm_r, nc, error, options, tolerance, has_not_converged, norm_q, w, itermax) \
                         shared(tmp_criteria1, tmp_criteria2, freeze_contacts, number_of_freezed_contact)
    {
    local_problem = (LinearComplementarityProblem *)malloc(sizeof(*local_problem));
    local_problem->M = NM_new();
    local_problem->M->storageType = NM_DENSE;
    local_problem->M->size0 = 2;
    local_problem->M->size1 = 2;
    local_problem->q = (double *)malloc(2 * sizeof(double));

    while ((iter < itermax) && has_not_converged) {
      // light_error_sum = 0.0;
      light_error_2 = 0.0;
      /* Loop over the rows of blocks in blmat */
      /* contact: current row (of blocks) number */
      
      #pragma omp single
      {
      number_of_freezed_contact = 0;
      tmp_criteria1 = tolerance * tolerance / (nc * nc * 10);
      tmp_criteria2 = *norm_r * *norm_r / (nc * nc * 1000);
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
        #pragma omp for reduction(+:light_error_sum)
        for (int v = 0; v < partition_size[color]; v++) {
          contact = partitions[color][v];

          if (freeze_contacts[contact] > 0) {
          /* we skip freeze contacts */
          freeze_contacts[contact] -= 1;
          continue;
          }

          pos = 2 * contact;
          localreaction[0] = z[pos];
          localreaction[1] = z[pos + 1];

          /* Local problem formalization */
          fc2d_nsgs_buildLocalProblem_parallel(contact, problem, diagonal_blocks, local_problem, z);

          /* Solve local problem */
          fc2d_nsgs_local_solve(local_problem->M->matrix0, diagonal_block_determinant[contact],
                                  local_problem->q, problem->mu[contact], localreaction);

          light_error_2 = light_error_squared(localreaction, &z[pos]);

          // #pragma omp atomic update
          light_error_sum += light_error_2;

          int relative_convergence_criteria =
              light_error_2 <= tmp_criteria1 * squared_norm(localreaction);
          int small_reaction_criteria = squared_norm(localreaction) <= tmp_criteria2;
          if ((relative_convergence_criteria || small_reaction_criteria) && iter >= 10) {
            /* we  freeze the contact for n iterations*/
            freeze_contacts[contact] = iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT];
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
        error = calculateLightError(light_error_sum, nc, z, norm_r);
        has_not_converged = determine_convergence(error, tolerance, iter, options);
      } else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                 SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
        error = calculateLightError(light_error_sum, nc, z, norm_r);
        has_not_converged = determine_convergence_with_full_final(
            problem, options, z, w, &tolerance, norm_q, error, iter);
      }
      light_error_sum = 0.;
      ++iter;
      }
    }  // end while loop
    free(local_problem->q);
    free(local_problem->M);
    free(local_problem);
  } // end parallel region 

  /***********************/
  /* NO FREEZING CONTACT */
  /***********************/
  } else {
    unsigned int contact;
    unsigned int pos;
    double light_error_sum = 0.;
    double localreaction[2];
    
    #pragma omp parallel default(none) \
                         private(contact, pos, local_problem, localreaction) \
                         shared(problem, diagonal_blocks, diagonal_block_determinant, z, partition_size, partitions, iter) \
                         shared(iparam, light_error_sum, n_colors, norm_r, nc, error, options, tolerance, has_not_converged, norm_q, w, itermax)
    {
    local_problem = (LinearComplementarityProblem *)malloc(sizeof(*local_problem));
    local_problem->M = NM_new();
    local_problem->M->storageType = NM_DENSE;
    local_problem->M->size0 = 2;
    local_problem->M->size1 = 2;
    local_problem->q = (double *)malloc(2 * sizeof(double));
    while ((iter < itermax) && has_not_converged) {
      // light_error_sum = 0.0;
      /* Loop over the rows of blocks in blmat */
      for (size_t color = 0; color < n_colors; color++) {
          #pragma omp for reduction(+:light_error_sum)
          for (int v = 0; v < partition_size[color]; v++) {
            contact = partitions[color][v];
            pos = 2 * contact;
            localreaction[0] = z[pos];
            localreaction[1] = z[pos + 1];

            /* Local problem formalization */
            fc2d_nsgs_buildLocalProblem_parallel(contact, problem, diagonal_blocks, local_problem, z);

            /* Solve local problem */
            fc2d_nsgs_local_solve(local_problem->M->matrix0, diagonal_block_determinant[contact],
                                    local_problem->q, problem->mu[contact], localreaction);

            light_error_sum += light_error_squared(localreaction, &z[pos]);
            // light_error_sum += (z[pos] - localreaction[0]) * (z[pos] - localreaction[0]) + (z[pos + 1] - localreaction[1]) * (z[pos + 1] - localreaction[1]);

            /* Why do we need the "if" here ? 
               If the condition is not true we don't even evaluate the error so it does not make sense... 
               So I guess the condition will always be met and therefore I should use a reduction here  can remove the "if"
               SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL seems to be the default
               
               I should use a reduction in the OMP pragma to aggregate light_error_sum without #pragma omp atomic update            
            */

            /* if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT ||
                iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
                    #pragma omp atomic update
                    light_error_sum += (z[pos] - localreaction[0]) * (z[pos] - localreaction[0]) + (z[pos + 1] - localreaction[1]) * (z[pos + 1] - localreaction[1]);
            } */
            z[pos] = localreaction[0];
            z[pos + 1] = localreaction[1];
          }
        }  // end for loop

      /*  error evaluation */
      #pragma omp single
      {
      if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
          SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT) {
        error = calculateLightError(light_error_sum, nc, z, norm_r);
        has_not_converged = determine_convergence(error, tolerance, iter, options);
      } else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                 SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
        error = calculateLightError(light_error_sum, nc, z, norm_r);
        has_not_converged = determine_convergence_with_full_final(
            problem, options, z, w, &tolerance, norm_q, error, iter);
      }
      light_error_sum = 0.;
      ++iter;
      }
    }  // end while loop
    free(local_problem->q);
    free(local_problem->M);
    free(local_problem);
    }
    }
  /* Full criterium */
  if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
      SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
    error = calculateFullErrorFinal(problem, options, z, w, tolerance, norm_q);

    has_not_converged =
        determine_convergence(error, dparam[SICONOS_DPARAM_TOL], iter, options);
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

  free(partition_size);
  for (size_t i = 0; i < n_colors; i++) free(partitions[i]);
  free(partitions);
}


/* *********** */
/* EXPERIMENTS */
/* *********** */

static void fc2d_nsgs_buildLocalProblem_2(unsigned int contact, FrictionContactProblem *problem,
                                        double *blocks_contiguous,
                                        double *diagonal_blocks,
					                              size_t *index1_data,
					                              size_t *index2_data,
                                        LinearComplementarityProblem *local_problem,
                                        double *reaction)
{
    local_problem->M->matrix0 = &diagonal_blocks[4 * contact];

    local_problem->M->size0 = 2; // Necessary ?
    local_problem->M->size1 = 2;

    local_problem->q[0] = problem->q[contact * 2];
    local_problem->q[1] = problem->q[contact * 2 + 1];

    size_t colNumber;
    // size_t nextColNumber;
    // unsigned int prefetch_distance = 4;

    for (size_t blockNum = index1_data[contact]; blockNum < index1_data[contact + 1]; blockNum++)
    {
      
      colNumber = index2_data[blockNum];
      // nextColNumber = index2_data[blockNum + prefetch_distance];
      // __builtin_prefetch(&reaction[2 * nextColNumber], 0, 3);

      if (colNumber != contact)
      {
        mvp2x2(&blocks_contiguous[blockNum * 4], &reaction[2 * colNumber], local_problem->q);
      }

    }
}

static double *fc2d_extract_diagonal_blocks(FrictionContactProblem *problem)
{
    unsigned int nc = problem->numberOfContacts;
    double *sbcm = (double *)calloc(4 * nc, sizeof(double));
    double *block;
    int diagPos;

    for (unsigned int i = 0; i < nc; ++i)
    {
        diagPos = SBM_diagonal_block_index(problem->M->matrix1, i);
        block = problem->M->matrix1->block[diagPos];
        for (int j = 0; j < 4; j++)
            sbcm[i * 4 + j] = block[j];
    }

    return sbcm;
}

static double *fc2d_nsgs_compute_local_problem_determinant_2(unsigned int nc, double *diagonal_blocks)
{
    double *diagonal_block_determinant = (double *)calloc(sizeof(double), nc);
    double *block;

    for (unsigned int i = 0; i < nc; ++i)
    {
        block = &diagonal_blocks[4 * i];
        diagonal_block_determinant[i] = block[0] * block[3] - block[1] * block[2];
        if (diagonal_block_determinant[i] < DBL_EPSILON)
        {
            free(diagonal_block_determinant);
            return NULL;
        }
    }
    return diagonal_block_determinant;
}

/* MAIN FUNCTION */
void test_solver(FrictionContactProblem *problem, double *z, double *w, int *info, SolverOptions *options)
{
    unsigned int nc = problem->numberOfContacts;
    unsigned int n = problem->numberOfContacts * problem->dimension;

    printf("Number of contacts = %d, dimension = %d\n", nc, n / nc);
    
    /* Coloring */
    size_t n_colors = 0;
    size_t *partition_size = NULL;
    size_t **partitions = NULL;
    size_t *inv_permutation = (size_t *)malloc((size_t)nc * sizeof(size_t));
    color_graph_block_permut(problem->numberOfContacts, problem->M, &n_colors, &partition_size, &partitions, inv_permutation);

    unsigned int *sum_sizes = (unsigned int *)calloc(n_colors + 1, sizeof(unsigned int));
    for (unsigned int color = 0; color < n_colors; color++)
    {
        sum_sizes[color + 1] = sum_sizes[color] + partition_size[color];
    }

    /* Permutate rows and columns of SBM */
    SparseBlockStructuredMatrix *SBM_col_permuted = SBM_new();
    SparseBlockStructuredMatrix *SBM_permuted = SBM_new();
    unsigned int *rowIndex = (unsigned int *)malloc((unsigned int)nc * sizeof(unsigned int));
    for (int i = 0; i < nc; i++)
        rowIndex[inv_permutation[i]] = i;

    SBM_column_permutation(rowIndex, problem->M->matrix1, SBM_col_permuted);
    SBM_row_permutation_copy(inv_permutation, SBM_col_permuted, SBM_permuted);

    problem->M->matrix1 = SBM_permuted;

    /* Cool stuff */
    /* unsigned int *columns_seen = (unsigned int *)calloc(nc, sizeof(unsigned int));
    size_t col;
    unsigned int nbdiffcol = 0, already_seen = 0;

    for (size_t color = 0; color < n_colors; color++)
    {

      printf("\n\nCOLOR %d\n\n", color);

      for (unsigned int contact = sum_sizes[color]; contact < sum_sizes[color + 1]; contact++)
      {

        for (size_t blockNum = SBM_permuted->index1_data[contact]; blockNum < SBM_permuted->index1_data[contact + 1]; blockNum++)
        {

          col = SBM_permuted->index2_data[blockNum];
          printf(" %d ", col);
          columns_seen[col]++;
          if (columns_seen[col] == 1)
          {
            nbdiffcol++;
          }
          else
          {
            already_seen++;
            // printf(" : SEEN ");
          }

        }

        printf("\n");
      }

      // printf("Color %d, number of different columns : %d / %d, already_seen = %d\n", color, nbdiffcol, nc, already_seen);

      for (unsigned int i = 0; i < nc; i++) columns_seen[i] = 0;
      nbdiffcol = 0;
      already_seen = 0;
    } */

    /* Store all blocks in contiguous array */
    unsigned int nbblocks = problem->M->matrix1->nbblocks;
    double *blocks_contiguous = (double *)malloc(nbblocks * 4 * sizeof(double));
    double *current_block;
    for (unsigned int blockNum = 0; blockNum < nbblocks; blockNum++)
    {
        current_block = problem->M->matrix1->block[blockNum];
        for (unsigned int j = 0; j < 4; j++)
        {
            blocks_contiguous[blockNum * 4 + j] = current_block[j];
        }
    }

    double *diagonal_blocks = fc2d_extract_diagonal_blocks(problem);
    double *diagonal_block_determinant = fc2d_nsgs_compute_local_problem_determinant_2(nc, diagonal_blocks);

    LinearComplementarityProblem *local_problem = (LinearComplementarityProblem *)malloc(sizeof(LinearComplementarityProblem));
    local_problem->M = NM_new();
    local_problem->M->storageType = NM_DENSE;
    local_problem->M->size0 = 2;
    local_problem->M->size1 = 2;
    local_problem->q = (double *)malloc(2 * sizeof(double));

    double localreaction[2];
    unsigned int pos;
    unsigned int contact;

    unsigned int permuted_contact;

    double *mu_permuted = (double *)malloc(nc * sizeof(double));
    double *q_permuted = (double *)malloc(nc * 2 * sizeof(double));
    for (unsigned int i = 0; i < nc; i++)
    {
        q_permuted[2 * i] = problem->q[2 * inv_permutation[i]];
        q_permuted[2 * i + 1] = problem->q[2 * inv_permutation[i] + 1];
        mu_permuted[i] = problem->mu[inv_permutation[i]];
    }

    problem->q = q_permuted;
    problem->mu = mu_permuted;

    double norm_q = cblas_dnrm2(nc * 2, problem->q, 1);
    double norm_r[] = {INFINITY};

    for (int i = 0; i < n; i++)
        z[i] = 1.;

    double error = INFINITY;
    double light_error_sum = 0.;
    int has_not_converged = 1;
    int iter = 0;

    int *iparam = options->iparam;
    double tolerance = 1e-3;
    int n_iter = iparam[SICONOS_IPARAM_MAX_ITER];

    size_t *index1_data;
    size_t *index2_data;

    unsigned int max_blocks_in_one_line = 0;

    for (unsigned int line = 0; line < problem->numberOfContacts; line++)
    {
      if (max_blocks_in_one_line < (problem->M->matrix1->index1_data[line + 1] - problem->M->matrix1->index1_data[line]))
      {
        max_blocks_in_one_line = (problem->M->matrix1->index1_data[line + 1] - problem->M->matrix1->index1_data[line]);
      }
    }
    printf("Maximum number of non empty blocks in one line = %d\n", max_blocks_in_one_line);

    double time = omp_get_wtime();

    #pragma omp parallel default(none) \
                         private(pos, local_problem, localreaction, index1_data, index2_data) \
                         shared(problem, diagonal_blocks, diagonal_block_determinant, blocks_contiguous, z, sum_sizes, n_colors, n_iter) \
                         // shared(light_error_sum, iparam, error, norm_r, norm_q, tolerance, has_not_converged, iter, options, w, nc)
    {

    printf("Thread num: %d\n", omp_get_thread_num());

    // Copy index1_data and index2_data
    index1_data = (size_t *)malloc((problem->numberOfContacts + 1) * sizeof(size_t));
    index2_data = (size_t *)malloc(problem->M->matrix1->nbblocks * sizeof(size_t));

    memcpy(index1_data, problem->M->matrix1->index1_data, (problem->numberOfContacts + 1) * sizeof(size_t));
    memcpy(index2_data, problem->M->matrix1->index2_data, problem->M->matrix1->nbblocks * sizeof(size_t));

    double *z_local = (double *)malloc(problem->numberOfContacts * 2 * sizeof(double));
    memcpy(z_local, z, problem->numberOfContacts * 2 * sizeof(double));

    double *blocks_contiguous_local = (double *)malloc(problem->M->matrix1->nbblocks * 4 * sizeof(double));
    memcpy(blocks_contiguous_local, blocks_contiguous, problem->M->matrix1->nbblocks * 4 * sizeof(double));
    double *diagonal_blocks_local = fc2d_extract_diagonal_blocks(problem);

    // Allocate local_problem for each thread
    local_problem = (LinearComplementarityProblem *)malloc(sizeof(*local_problem));
    local_problem->M = NM_new();
    local_problem->M->storageType = NM_DENSE;
    local_problem->M->size0 = 2;
    local_problem->M->size1 = 2;
    local_problem->q = (double *)malloc(2 * sizeof(double));
  
    for (int i = 0; i < n_iter; i++)
    {
        for (size_t color = 0; color < n_colors; color++)
        {
	    // #pragma omp for reduction(+:light_error_sum)
            #pragma omp for
            for (unsigned int permuted_contact = sum_sizes[color]; permuted_contact < sum_sizes[color + 1]; permuted_contact++)
            {
                // printf("[%d] Contact : %d\n", omp_get_thread_num(), permuted_contact);
                // contact = partitions[color][v];
                pos = 2 * permuted_contact;
                localreaction[0] = z[pos];
                localreaction[1] = z[pos + 1];

                fc2d_nsgs_buildLocalProblem_2(permuted_contact, problem, blocks_contiguous_local, diagonal_blocks_local, index1_data, index2_data, local_problem, z_local);

                /* printf("Contact %d, q[0] = %e, q[1] = %e\n", contact, local_problem->q[0], local_problem->q[1]);
                printf("M->matrix0 = %e %e %e %e \n", local_problem->M->matrix0[0], local_problem->M->matrix0[1], local_problem->M->matrix0[2], local_problem->M->matrix0[3]);
                printf("%e %e\n", diagonal_block_determinant[permuted_contact], problem->mu[permuted_contact]);
                printf("%e %e\n", localreaction[0], localreaction[1]); */

                /* Solve local problem */
                // fc2d_nsgs_local_solve(local_problem->M->matrix0, diagonal_block_determinant[permuted_contact], local_problem->q, problem->mu[permuted_contact], localreaction);

                // light_error_sum += light_error_squared(localreaction, &z[pos]);

                /* printf("AFTER LOCAL SOLVE\n");
                printf("%e %e\n", localreaction[0], localreaction[1]); */
                // z[pos] = localreaction[0];
                // z[pos + 1] = localreaction[1];
                // printf("\n");
            }

            // printf("\nCOLOR %d DONE\n", color);
        }

	/* 
	#pragma omp single
        {
        if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT) {
            error = calculateLightError(light_error_sum, nc, z, norm_r);
            has_not_converged = determine_convergence(error, tolerance, iter, options);
        } else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] == SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
            error = calculateLightError(light_error_sum, nc, z, norm_r);
            has_not_converged = determine_convergence_with_full_final(problem, options, z, w, &tolerance, norm_q, error, iter);
        }

	light_error_sum = 0.0;
        iter++;
	}
	*/

        /* for (int j = 0; j < n; j++)
            printf("%e \n", z[j]); */
    }
    // End of parallel section
    }

    /* for (int j = 0; j < n; j++)
        printf("%e \n", z[j]); */

    time = omp_get_wtime() - time;
    printf("Time iterating: %e seconds\n", time);

    double *z_permut = (double *)malloc(2 * nc * sizeof(double));
    for (unsigned int i = 0; i < nc; i++)
    {
        z_permut[2 * inv_permutation[i]] = z[2 * i];
        z_permut[2 * inv_permutation[i] + 1] = z[2 * i + 1];
    }

    //for (int j = 0; j < n; j++)
    //    printf("%e \n", z_permut[j]);
}