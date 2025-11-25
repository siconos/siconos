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

void fc2d_sbm_to_bsr(SparseBlockStructuredMatrix *sbm, fc2d_cuda_context_t *ctx) {
  unsigned int nc = sbm->blocknumber0;
  unsigned int nbblocks = sbm->nbblocks;

  /* Allocate host arrays for BSR format */
  int *h_bsrRowPtr = (int *)malloc((nc + 1) * sizeof(int));
  int *h_bsrColInd = (int *)malloc(nbblocks * sizeof(int));
  double *h_bsrVal = (double *)malloc(nbblocks * 4 * sizeof(double));

  if(!h_bsrRowPtr || !h_bsrColInd || !h_bsrVal) {
    fprintf(stderr, "Failed to allocate host memory for BSR conversion\n");
    exit(EXIT_FAILURE);
  }

  /* Convert index1_data to bsrRowPtr (size_t to int) */
  for(size_t i = 0; i <= nc; i++) {
    h_bsrRowPtr[i] = (int)sbm->index1_data[i];
  }

  /* Convert index2_data to bsrColInd and blocks to bsrVal */
  for(size_t i = 0; i < nbblocks; i++) {
    h_bsrColInd[i] = (int)sbm->index2_data[i];

    /* Copy 2x2 block (already in column-major format) */
    double *block = sbm->block[i];
    h_bsrVal[i * 4 + 0] = block[0]; /* column 0, row 0 */
    h_bsrVal[i * 4 + 1] = block[1]; /* column 0, row 1 */
    h_bsrVal[i * 4 + 2] = block[2]; /* column 1, row 0 */
    h_bsrVal[i * 4 + 3] = block[3]; /* column 1, row 1 */
  }

  /* Allocate device memory */
  CUDA_CHECK(cudaMalloc(&ctx->d_bsrRowPtr, (nc + 1) * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&ctx->d_bsrColInd, nbblocks * sizeof(int)));
  CUDA_CHECK(cudaMalloc(&ctx->d_bsrVal, nbblocks * 4 * sizeof(double)));

  /* Copy to device */
  CUDA_CHECK(cudaMemcpy(ctx->d_bsrRowPtr, h_bsrRowPtr, (nc + 1) * sizeof(int),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(
      cudaMemcpy(ctx->d_bsrColInd, h_bsrColInd, nbblocks * sizeof(int), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(ctx->d_bsrVal, h_bsrVal, nbblocks * 4 * sizeof(double),
                        cudaMemcpyHostToDevice));

  /* Create cuSparse BSR matrix descriptor */
  CUSPARSE_CHECK(cusparseCreateBsr(&ctx->bsr_matrix,
                                   nc,                       /* rows (in blocks) */
                                   nc,                       /* cols (in blocks) */
                                   nbblocks,                 /* nnz (number of blocks) */
                                   ctx->d_bsrRowPtr,         /* bsrRowPtr */
                                   ctx->d_bsrColInd,         /* bsrColInd */
                                   ctx->d_bsrVal,            /* bsrVal */
                                   2,                        /* block dimension */
                                   CUSPARSE_INDEX_32I,       /* row pointer type */
                                   CUSPARSE_INDEX_32I,       /* col index type */
                                   CUSPARSE_INDEX_BASE_ZERO, /* base index */
                                   CUDA_R_64F                /* data type */
                                   ));

  /* Free host memory */
  free(h_bsrRowPtr);
  free(h_bsrColInd);
  free(h_bsrVal);
}
