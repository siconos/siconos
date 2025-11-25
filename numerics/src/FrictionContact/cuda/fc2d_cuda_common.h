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

#ifndef FC2D_CUDA_COMMON_H
#define FC2D_CUDA_COMMON_H

#include <cuda_runtime.h>
#include <cusparse.h>
#include <stdio.h>
#include <stdlib.h>

#include "FrictionContactProblem.h"
#include "SolverOptions.h"
#include "SparseBlockMatrix.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/**
 * GPU context for fc2d CUDA solvers.
 * Contains all GPU-allocated data and cuSparse structures.
 */
typedef struct {
  /* cuSparse handle */
  cusparseHandle_t cusparse_handle;

  /* Matrix in BSR format */
  cusparseSpMatDescr_t bsr_matrix;
  int *d_bsrRowPtr;   /* Device: row pointers (size nc+1) */
  int *d_bsrColInd;   /* Device: column indices (size nbblocks) */
  double *d_bsrVal;   /* Device: block values (size nbblocks*4) */

  /* Diagonal blocks and determinants */
  double *d_diagonal_blocks;       /* Device: 2x2 diagonal blocks (size 4*nc) */
  double *d_diagonal_determinants; /* Device: determinants (size nc) */

  /* Problem vectors */
  double *d_z;  /* Device: reaction vector (size 2*nc) */
  double *d_q;  /* Device: right-hand side (size 2*nc) */
  double *d_mu; /* Device: friction coefficients (size nc) */

  /* Graph coloring data */
  size_t n_colors;
  size_t *d_sum_sizes; /* Device: color boundaries (size n_colors+1) */

  /* Workspace for error computation */
  double *d_error_local; /* Device: per-thread error contributions */
  double *d_norm_local;  /* Device: per-thread norm contributions */

  /* Problem size */
  unsigned int nc;       /* Number of contacts */
  unsigned int nbblocks; /* Number of non-null blocks */

} fc2d_cuda_context_t;

/**
 * Error checking macros for CUDA and cuSparse calls
 */
#define CUDA_CHECK(call)                                                                         \
  do {                                                                                           \
    cudaError_t err = call;                                                                      \
    if(err != cudaSuccess) {                                                                     \
      fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,                          \
              cudaGetErrorString(err));                                                          \
      exit(EXIT_FAILURE);                                                                        \
    }                                                                                            \
  } while(0)

#define CUSPARSE_CHECK(call)                                                                     \
  do {                                                                                           \
    cusparseStatus_t err = call;                                                                 \
    if(err != CUSPARSE_STATUS_SUCCESS) {                                                         \
      fprintf(stderr, "cuSPARSE error at %s:%d: %d\n", __FILE__, __LINE__, err);                \
      exit(EXIT_FAILURE);                                                                        \
    }                                                                                            \
  } while(0)

/**
 * Context management functions
 */

/**
 * Create and initialize a CUDA context for fc2d solvers.
 * Allocates all necessary GPU memory.
 *
 * \param nc Number of contacts
 * \param nbblocks Number of non-null blocks in the matrix
 * \return Pointer to initialized context, or NULL on error
 */
fc2d_cuda_context_t *fc2d_cuda_context_create(unsigned int nc, unsigned int nbblocks);

/**
 * Destroy a CUDA context and free all associated GPU memory.
 *
 * \param ctx Context to destroy
 */
void fc2d_cuda_context_destroy(fc2d_cuda_context_t *ctx);

/**
 * Matrix conversion functions
 */

/**
 * Convert a SparseBlockStructuredMatrix to BSR format and transfer to GPU.
 * Creates cuSparse BSR matrix descriptor.
 *
 * \param sbm Source sparse block matrix
 * \param ctx CUDA context (will be populated with BSR data)
 */
void fc2d_sbm_to_bsr(SparseBlockStructuredMatrix *sbm, fc2d_cuda_context_t *ctx);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif /* FC2D_CUDA_COMMON_H */
