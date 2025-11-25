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

#include "fc2d_Solvers.h"
#include "fc2d_cuda_common.h"
#include "fc2d_cuda_kernels.h"

#include "Friction_cst.h"
#include "SparseBlockMatrix.h"
#include "graph_tools.h"
#include "numerics_verbose.h"

#include <float.h>
#include <math.h>
#include <string.h>

/**
 * GPU-accelerated NSGS solver with graph coloring and permutation.
 */
void fc2d_nsgs_graph_permut_cuda(FrictionContactProblem *problem, double *z, double *w, int *info,
                                  SolverOptions *options) {
  /* Get solver parameters */
  unsigned int nc = problem->numberOfContacts;
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = options->dparam[SICONOS_DPARAM_TOL];

  /* STEP 1: Graph coloring and permutation (on CPU, same as original) */
  size_t n_colors = 0;
  size_t *sum_sizes = NULL;
  size_t *inv_permutation = (size_t *)malloc(nc * sizeof(size_t));

  color_graph_block_permut(nc, problem->M, &n_colors, &sum_sizes, inv_permutation);

  numerics_printf("-- FC2D NSGS CUDA - Graph colored with %zu colors\n", n_colors);

  /* STEP 2: Permute matrix */
  SparseBlockStructuredMatrix *SBM_col_permuted = SBM_new();
  SparseBlockStructuredMatrix *SBM_permuted = SBM_new();
  unsigned int *rowIndex = (unsigned int *)malloc(nc * sizeof(unsigned int));

  for(unsigned int i = 0; i < nc; i++)
    rowIndex[inv_permutation[i]] = i;

  SBM_column_permutation(rowIndex, problem->M->matrix1, SBM_col_permuted);
  SBM_row_permutation_copy(inv_permutation, SBM_col_permuted, SBM_permuted);
  free(rowIndex);

  /* Temporarily swap matrix */
  SparseBlockStructuredMatrix *old_matrix1 = problem->M->matrix1;
  problem->M->matrix1 = SBM_permuted;

  /* STEP 3: Extract and permute problem data */
  unsigned int nbblocks = SBM_permuted->nbblocks;

  /* Extract diagonal blocks and compute determinants */
  double *h_diagonal_blocks = (double *)malloc(4 * nc * sizeof(double));
  double *h_diagonal_determinants = (double *)malloc(nc * sizeof(double));

  for(unsigned int i = 0; i < nc; i++) {
    int diagPos = SBM_diagonal_block_index(SBM_permuted, i);
    double *block = SBM_permuted->block[diagPos];

    /* Copy block */
    for(int j = 0; j < 4; j++)
      h_diagonal_blocks[i * 4 + j] = block[j];

    /* Compute determinant */
    h_diagonal_determinants[i] = block[0] * block[3] - block[1] * block[2];

    if(h_diagonal_determinants[i] < DBL_EPSILON) {
      fprintf(stderr, "FC2D NSGS CUDA: Singular diagonal block at contact %u\n", i);
      *info = 2;
      goto cleanup;
    }
  }

  /* Permute q, mu, and z */
  double *q_permuted = (double *)malloc(2 * nc * sizeof(double));
  double *mu_permuted = (double *)malloc(nc * sizeof(double));
  double *z_permuted = (double *)malloc(2 * nc * sizeof(double));

  for(unsigned int i = 0; i < nc; i++) {
    q_permuted[2 * i] = problem->q[2 * inv_permutation[i]];
    q_permuted[2 * i + 1] = problem->q[2 * inv_permutation[i] + 1];
    mu_permuted[i] = problem->mu[inv_permutation[i]];
    z_permuted[2 * i] = z[2 * inv_permutation[i]];
    z_permuted[2 * i + 1] = z[2 * inv_permutation[i] + 1];
  }

  /* STEP 4: Create GPU context and transfer data */
  fc2d_cuda_context_t *ctx = fc2d_cuda_context_create(nc, nbblocks);
  if(!ctx) {
    fprintf(stderr, "FC2D NSGS CUDA: Failed to create GPU context\n");
    *info = 2;
    goto cleanup;
  }

  ctx->n_colors = n_colors;

  /* Convert SBM to BSR and transfer to GPU */
  fc2d_sbm_to_bsr(SBM_permuted, ctx);

  /* Transfer vectors to GPU */
  CUDA_CHECK(cudaMemcpy(ctx->d_diagonal_blocks, h_diagonal_blocks, 4 * nc * sizeof(double),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(ctx->d_diagonal_determinants, h_diagonal_determinants, nc * sizeof(double),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(ctx->d_q, q_permuted, 2 * nc * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(ctx->d_mu, mu_permuted, nc * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(ctx->d_z, z_permuted, 2 * nc * sizeof(double), cudaMemcpyHostToDevice));

  /* Transfer color boundaries */
  CUDA_CHECK(cudaMalloc(&ctx->d_sum_sizes, (n_colors + 1) * sizeof(size_t)));
  CUDA_CHECK(cudaMemcpy(ctx->d_sum_sizes, sum_sizes, (n_colors + 1) * sizeof(size_t),
                        cudaMemcpyHostToDevice));

  numerics_printf("-- FC2D NSGS CUDA - Data transferred to GPU\n");

  /* STEP 5: Iterative solver loop */
  int iter = 0;
  double error = INFINITY;
  int has_not_converged = 1;

  /* Copy sum_sizes to host for iteration */
  size_t *h_sum_sizes = (size_t *)malloc((n_colors + 1) * sizeof(size_t));
  memcpy(h_sum_sizes, sum_sizes, (n_colors + 1) * sizeof(size_t));

  numerics_printf("-- FC2D NSGS CUDA - Starting iterations (max: %d, tol: %e)\n", itermax,
                  tolerance);

  while(iter < itermax && has_not_converged) {
    double error_sum = 0.0;
    double norm_sum = 0.0;

    /* Loop over colors (sequential) */
    for(size_t color = 0; color < n_colors; color++) {
      unsigned int contact_start = (unsigned int)h_sum_sizes[color];
      unsigned int contact_end = (unsigned int)h_sum_sizes[color + 1];
      unsigned int num_contacts = contact_end - contact_start;

      if(num_contacts == 0)
        continue;

      /* Launch kernel for this color */
      fc2d_nsgs_color_iteration(contact_start, contact_end, ctx->d_bsrRowPtr, ctx->d_bsrColInd,
                                ctx->d_bsrVal, ctx->d_diagonal_blocks,
                                ctx->d_diagonal_determinants, ctx->d_q, ctx->d_mu, ctx->d_z,
                                ctx->d_error_local, ctx->d_norm_local,
                                0 /* default stream */);

      /* Reduce error and norm for this color */
      double color_error, color_norm;
      reduce_error_norm(ctx->d_error_local, ctx->d_norm_local, num_contacts, &color_error,
                        &color_norm);

      error_sum += color_error;
      norm_sum += color_norm;
    }

    /* Compute error */
    error = sqrt(error_sum);
    double norm_r = sqrt(norm_sum);

    if(fabs(norm_r) > DBL_EPSILON) {
      error /= norm_r;
    }

    /* Check convergence */
    if(error < tolerance) {
      has_not_converged = 0;
      numerics_printf("-- FC2D NSGS CUDA - Iteration %i Residual = %14.7e < %7.3e (converged)\n",
                      iter, error, tolerance);
    } else {
      if(verbose > 0 && (iter % 10 == 0 || iter < 5)) {
        numerics_printf("-- FC2D NSGS CUDA - Iteration %i Residual = %14.7e > %7.3e\n", iter,
                        error, tolerance);
      }
    }

    iter++;
  }

  numerics_printf("-- FC2D NSGS CUDA - Finished after %d iterations (error = %e)\n", iter, error);

  /* STEP 6: Copy results back to host */
  CUDA_CHECK(cudaMemcpy(z_permuted, ctx->d_z, 2 * nc * sizeof(double), cudaMemcpyDeviceToHost));

  /* STEP 7: Unpermute results */
  for(unsigned int i = 0; i < nc; i++) {
    z[2 * inv_permutation[i]] = z_permuted[2 * i];
    z[2 * inv_permutation[i] + 1] = z_permuted[2 * i + 1];
  }

  /* Set output parameters */
  *info = has_not_converged;
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  options->dparam[SICONOS_DPARAM_RESIDU] = error;

cleanup:
  /* Cleanup GPU context */
  if(ctx)
    fc2d_cuda_context_destroy(ctx);

  /* Free host memory */
  if(h_diagonal_blocks)
    free(h_diagonal_blocks);
  if(h_diagonal_determinants)
    free(h_diagonal_determinants);
  if(q_permuted)
    free(q_permuted);
  if(mu_permuted)
    free(mu_permuted);
  if(z_permuted)
    free(z_permuted);
  if(h_sum_sizes)
    free(h_sum_sizes);
  if(inv_permutation)
    free(inv_permutation);
  if(sum_sizes)
    free(sum_sizes);

  /* Restore original matrix */
  problem->M->matrix1 = old_matrix1;
  if(SBM_col_permuted) {
    SBM_clear(SBM_col_permuted);
    free(SBM_col_permuted);
  }
  if(SBM_permuted) {
    SBM_clear(SBM_permuted);
    free(SBM_permuted);
  }
}
