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

#include "fc2d_cuda_common.h"
#include <string.h>

fc2d_cuda_context_t *fc2d_cuda_context_create(unsigned int nc, unsigned int nbblocks) {
  fc2d_cuda_context_t *ctx = (fc2d_cuda_context_t *)malloc(sizeof(fc2d_cuda_context_t));

  if(!ctx) {
    fprintf(stderr, "Failed to allocate CUDA context\n");
    return NULL;
  }

  /* Initialize structure */
  memset(ctx, 0, sizeof(fc2d_cuda_context_t));
  ctx->nc = nc;
  ctx->nbblocks = nbblocks;

  /* Create cuSPARSE handle */
  CUSPARSE_CHECK(cusparseCreate(&ctx->cusparse_handle));

  /* Allocate device memory for vectors */
  CUDA_CHECK(cudaMalloc(&ctx->d_z, 2 * nc * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&ctx->d_q, 2 * nc * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&ctx->d_mu, nc * sizeof(double)));

  /* Allocate device memory for diagonal blocks and determinants */
  CUDA_CHECK(cudaMalloc(&ctx->d_diagonal_blocks, 4 * nc * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&ctx->d_diagonal_determinants, nc * sizeof(double)));

  /* Allocate workspace for error computation */
  CUDA_CHECK(cudaMalloc(&ctx->d_error_local, nc * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&ctx->d_norm_local, nc * sizeof(double)));

  /* BSR matrix pointers will be allocated in fc2d_sbm_to_bsr() */
  ctx->d_bsrRowPtr = NULL;
  ctx->d_bsrColInd = NULL;
  ctx->d_bsrVal = NULL;
  ctx->bsr_matrix = NULL;

  /* Color boundaries will be set by solver */
  ctx->d_sum_sizes = NULL;
  ctx->n_colors = 0;

  return ctx;
}

void fc2d_cuda_context_destroy(fc2d_cuda_context_t *ctx) {
  if(!ctx)
    return;

  /* Destroy cuSPARSE objects */
  if(ctx->bsr_matrix) {
    cusparseDestroySpMat(ctx->bsr_matrix);
  }
  cusparseDestroy(ctx->cusparse_handle);

  /* Free device memory */
  if(ctx->d_bsrRowPtr)
    cudaFree(ctx->d_bsrRowPtr);
  if(ctx->d_bsrColInd)
    cudaFree(ctx->d_bsrColInd);
  if(ctx->d_bsrVal)
    cudaFree(ctx->d_bsrVal);
  if(ctx->d_z)
    cudaFree(ctx->d_z);
  if(ctx->d_q)
    cudaFree(ctx->d_q);
  if(ctx->d_mu)
    cudaFree(ctx->d_mu);
  if(ctx->d_diagonal_blocks)
    cudaFree(ctx->d_diagonal_blocks);
  if(ctx->d_diagonal_determinants)
    cudaFree(ctx->d_diagonal_determinants);
  if(ctx->d_sum_sizes)
    cudaFree(ctx->d_sum_sizes);
  if(ctx->d_error_local)
    cudaFree(ctx->d_error_local);
  if(ctx->d_norm_local)
    cudaFree(ctx->d_norm_local);

  /* Free context structure */
  free(ctx);
}
