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
#include <limits.h>
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

/* CUDA */
#include <cuda_runtime.h>
#include <cusparse.h>

#define CHECK_CUDA(func)                                                   \
  {                                                                        \
    cudaError_t status = (func);                                           \
    if (status != cudaSuccess) {                                           \
      printf("CUDA API failed at line %d with error: %s (%d)\n", __LINE__, \
             cudaGetErrorString(status), status);                          \
      *info = 1;                                                           \
      return;                                                              \
    }                                                                      \
  }

#define CHECK_CUSPARSE(func)                                                   \
  {                                                                            \
    cusparseStatus_t status = (func);                                          \
    if (status != CUSPARSE_STATUS_SUCCESS) {                                   \
      printf("CUSPARSE API failed at line %d with error: %s (%d)\n", __LINE__, \
             cusparseGetErrorString(status), status);                          \
      *info = 1;                                                               \
      return;                                                                  \
    }                                                                          \
  }

#define SGN(x) ((x) < 0 ? -1 : (x) > 0 ? 1 : 0)

/* Extract diagonal blocks from a 2D problem, as a contiguous array.
 * Only works on SBM matrices.
 */
static double* fc2d_extract_diagonal_blocks(FrictionContactProblem* problem) {
  unsigned int nc = problem->numberOfContacts;
  double* sbcm = (double*)calloc(4 * nc, sizeof(double));
  double* block;
  int diagPos;

  for (unsigned int i = 0; i < nc; ++i) {
    diagPos = SBM_diagonal_block_index(problem->M->matrix1, i);
    block = problem->M->matrix1->block[diagPos];
    for (int j = 0; j < 4; j++) sbcm[i * 4 + j] = block[j];
  }

  return sbcm;
}

/* Free diagonal blocks array */
static void* fc2d_free_diagonal_blocks(double* sbcm) {
  free(sbcm);
  return NULL;
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
                                                 double* device_reaction, double* velocity,
                                                 double* tolerance, double norm_q,
                                                 double error, unsigned int iter) {
  int has_not_converged = 1;
  if (error < *tolerance) {
    // Copy z back to host
    cudaMemcpy(reaction, device_reaction, (2 * problem->numberOfContacts) * sizeof(double),
               cudaMemcpyDeviceToHost);
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

/* CUDA kernel to solve local problems,
 * compute and aggregate the errors.
 */
__global__ void fc2d_nsgs_local_solve_kernel_range_reduce_2(
    const double* W,        // Matrix blocks
    const double* D,        // Determinants
    const double* q_ax,     // Result from cuSPARSE (Ax)
    const double* q_const,  // The constant vector B
    const double* mu,       // Friction coefficients
    double* P,              // The solution vector z
    double* sumP2,          // Reduction: squared norm
    double* sumErr2,        // Reduction: squared error
    const int k1, const int k2) {
  extern __shared__ double shmem[];
  double* shP2 = shmem;
  double* shErr2 = shmem + blockDim.x;

  int tid = threadIdx.x;
  int local = blockIdx.x * blockDim.x + tid;
  int i = k1 + local;

  double localP2 = 0.0;
  double localErr2 = 0.0;

  if (i < k2) {
    // 1. FUSED READ: Compute q = Ax + B in registers
    // This replaces the cudaMemcpy(D2D)
    double qi0 = q_ax[2 * i] + q_const[2 * i];
    double qi1 = q_ax[2 * i + 1] + q_const[2 * i + 1];

    const double Di = D[i];
    const double mui = mu[i];
    const double* Wi = &W[4 * i];
    double* Pi = &P[2 * i];

    // Save old P for error calculation
    double Pold0 = Pi[0];
    double Pold1 = Pi[1];

    const double Wnn = Wi[0];
    const double Wtn = Wi[1];
    const double Wnt = Wi[2];
    const double Wtt = Wi[3];

    // 2. SOLVER LOGIC: Using our register-cached qi0 and qi1
    if (qi0 > 0.0) {
      Pi[0] = 0.0;
      Pi[1] = 0.0;
    } else {
      // Using the fused values here
      Pi[0] = -(Wtt * qi0 - Wnt * qi1) / Di;
      Pi[1] = -(-Wtn * qi0 + Wnn * qi1) / Di;

      double muPn = mui * Pi[0];

      if (fabs(Pi[1]) > muPn) {
        if (Pi[1] + muPn < 0.0) {
          Pi[0] = -qi0 / (Wnn - mui * Wnt);
          Pi[1] = -mui * Pi[0];
        } else {
          Pi[0] = -qi0 / (Wnn + mui * Wnt);
          Pi[1] = mui * Pi[0];
        }
      }
    }

    // 3. ACCUMULATION
    double dP0 = Pi[0] - Pold0;
    double dP1 = Pi[1] - Pold1;

    localP2 = Pi[0] * Pi[0] + Pi[1] * Pi[1];
    localErr2 = dP0 * dP0 + dP1 * dP1;
  }

  // --- Reduction Logic remains the same ---
  shP2[tid] = localP2;
  shErr2[tid] = localErr2;
  __syncthreads();

  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      shP2[tid] += shP2[tid + s];
      shErr2[tid] += shErr2[tid + s];
    }
    __syncthreads();
  }

  if (tid == 0) {
    atomicAdd(sumP2, shP2[0]);
    atomicAdd(sumErr2, shErr2[0]);
  }
}

void fc2d_nsgs_graph_permut_cuda_blocklegacy(FrictionContactProblem* problem, double* z,
                                             double* w, int* info, SolverOptions* options) {
  // Get solver parameters
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  unsigned int nc = problem->numberOfContacts;
  double norm_q = cblas_dnrm2(nc * 2, problem->q, 1);
  double norm_r = 0.0;
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = options->dparam[SICONOS_DPARAM_TOL];

  // Initialize cusparse
  cusparseHandle_t handle = NULL;
  CHECK_CUSPARSE(cusparseCreate(&handle))

  // Start coloring and permutation
  size_t n_colors = 0;
  size_t* sum_sizes = NULL;
  size_t* inv_permutation = (size_t*)malloc(nc * sizeof(size_t));

  color_graph_block_permut(nc, problem->M, &n_colors, &sum_sizes, inv_permutation);

  SparseBlockStructuredMatrix* SBM_problem = NULL;

  /* Create SBM if it does not exist */
  if (!problem->M->matrix1) {
    SBM_problem = SBM_new();
    switch (problem->M->storageType) {
      // Convert DENSE -> SBM
      case NM_DENSE: {
        SBM_from_dense(problem->dimension, problem->M->size1, problem->M->size0,
                       problem->M->matrix0, SBM_problem);
        break;
      }
      // Convert SPARSE -> SBM
      case NM_SPARSE: {
        // Get Csparse matrix
        CSparseMatrix* sparse;
        if (problem->M->matrix2->origin == NSM_CSR) {
          sparse = NM_csr(problem->M);
        } else {
          sparse = NM_csc_trans(problem->M);
        }

        SBM_from_csparse_2(problem->dimension, sparse, SBM_problem);
        break;
      }
    }
    problem->M->matrix1 = SBM_problem;
  }

  /* Switch storageType to SBM temporarily */
  NM_types old_storageType = problem->M->storageType;
  problem->M->storageType = NM_SPARSE_BLOCK;

  SparseBlockStructuredMatrix* SBM_col_permuted = SBM_new();
  SparseBlockStructuredMatrix* SBM_permuted = SBM_new();
  size_t* rowIndex = (size_t*)malloc(nc * sizeof(size_t));

  for (unsigned int i = 0; i < nc; i++) rowIndex[inv_permutation[i]] = i;

  SBM_column_permutation(rowIndex, problem->M->matrix1, SBM_col_permuted);
  SBM_row_permutation_copy(inv_permutation, SBM_col_permuted, SBM_permuted);
  free(rowIndex);

  SparseBlockStructuredMatrix* old_matrix1 = problem->M->matrix1;
  problem->M->matrix1 = SBM_permuted;

  /* Get diagonal blocks and determinants */
  double* diagonal_blocks = fc2d_extract_diagonal_blocks(problem);
  double* diagonal_block_determinant =
      fc2d_nsgs_compute_local_problem_determinant(nc, diagonal_blocks);

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
  // Permutation done

  // Create an array storing, for each color, the row offsets of the corresponding CSR
  // submatrix
  int* h_all_RowOffsets = (int*)malloc((nc + n_colors) * sizeof(int));
  size_t k = 0;
  for (unsigned int color = 0; color < n_colors; color++) {
    size_t start_line = sum_sizes[color];
    size_t end_line = sum_sizes[color + 1];

    for (unsigned int row = start_line; row < end_line + 1; row++) {
      h_all_RowOffsets[k] =
          SBM_permuted->index1_data[row] - SBM_permuted->index1_data[start_line];
      k++;
    }
  }

  // Copy this array to device
  int* d_all_RowOffsets = NULL;
  CHECK_CUDA(cudaMalloc(&d_all_RowOffsets, (nc + n_colors) * sizeof(int)))
  CHECK_CUDA(cudaMemcpy(d_all_RowOffsets, h_all_RowOffsets, (nc + n_colors) * sizeof(int),
                        cudaMemcpyHostToDevice))
  free(h_all_RowOffsets);

  // Column indices
  size_t nbblocks = SBM_permuted->nbblocks;
  int* h_all_ColIndices = (int*)malloc(nbblocks * sizeof(int));
  for (size_t b = 0; b < nbblocks; b++) h_all_ColIndices[b] = SBM_permuted->index2_data[b];
  int* d_all_ColIndices = NULL;

  CHECK_CUDA(cudaMalloc(&d_all_ColIndices,
                        nbblocks * sizeof(int)))  // column index for each block

  CHECK_CUDA(cudaMemcpy(d_all_ColIndices, h_all_ColIndices, nbblocks * sizeof(int),
                        cudaMemcpyHostToDevice))

  // Values in blocks
  double* d_all_Values = NULL;
  double* h_all_Values = (double*)malloc(
      4 * nbblocks * sizeof(double));  // to store all elements in a contiguous array
  double* current_block;
  int ii = 0;
  int diagPos = SBM_diagonal_block_index(SBM_permuted, ii);

  for (unsigned int blockNum = 0; blockNum < nbblocks; blockNum++) {
    if (blockNum == diagPos) {
      for (unsigned int j = 0; j < 4; j++) {
        h_all_Values[blockNum * 4 + j] = 0;
      }
      ii++;
      diagPos = SBM_diagonal_block_index(SBM_permuted, ii);
    } else {
      current_block = SBM_permuted->block[blockNum];
      for (unsigned int j = 0; j < 4; j++) {
        h_all_Values[blockNum * 4 + j] = current_block[j];
      }
    }
  }

  CHECK_CUDA(
      cudaMalloc(&d_all_Values, 4 * nbblocks * sizeof(double)))  // each block has 4 elements

  CHECK_CUDA(cudaMemcpy(d_all_Values, h_all_Values, 4 * nbblocks * sizeof(double),
                        cudaMemcpyHostToDevice))

  free(h_all_Values);

  // Create matrix descriptor
  cusparseMatDescr_t descr = 0;
  cusparseCreateMatDescr(&descr);
  cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
  cusparseSetMatType(descr,
                     CUSPARSE_MATRIX_TYPE_GENERAL);  // matrix type: general, symmmetric,
                                                     // Hermitian, triangular

  // Create vector for z, the reaction, which is the input of the product
  double* d_z = NULL;
  CHECK_CUDA(cudaMalloc(&d_z, (2 * nc) * sizeof(double)))
  CHECK_CUDA(cudaMemcpy(d_z, z, (2 * nc) * sizeof(double), cudaMemcpyHostToDevice))

  // Create the res vectors, which are the q's of the local problems
  double* d_q_new = NULL;
  CHECK_CUDA(cudaMalloc(&d_q_new, (2 * nc) * sizeof(double)))
  CHECK_CUDA(cudaMemset(d_q_new, 0., (2 * nc) * sizeof(double)))

  // Copy diagonal blocks to device
  double* d_diagonal_blocks = NULL;
  CHECK_CUDA(cudaMalloc(&d_diagonal_blocks, (4 * nc) * sizeof(double)))
  CHECK_CUDA(cudaMemcpy(d_diagonal_blocks, diagonal_blocks, (4 * nc) * sizeof(double),
                        cudaMemcpyHostToDevice))

  // We don't need diagonal blocks anymore
  fc2d_free_diagonal_blocks(diagonal_blocks);

  // Copy determinants of diagonal blocks to device
  double* d_determinants = NULL;
  CHECK_CUDA(cudaMalloc(&d_determinants, nc * sizeof(double)))
  CHECK_CUDA(cudaMemcpy(d_determinants, diagonal_block_determinant, nc * sizeof(double),
                        cudaMemcpyHostToDevice))

  // We don't need diagonal block determinants anymore
  free(diagonal_block_determinant);

  // Copy mu to device
  double* d_mu = NULL;
  CHECK_CUDA(cudaMalloc(&d_mu, nc * sizeof(double)))
  CHECK_CUDA(cudaMemcpy(d_mu, problem->mu, nc * sizeof(double), cudaMemcpyHostToDevice))

  // Constant q on device
  double* d_q_const = NULL;
  CHECK_CUDA(cudaMalloc(&d_q_const, (2 * nc) * sizeof(double)))
  CHECK_CUDA(
      cudaMemcpy(d_q_const, problem->q, (2 * nc) * sizeof(double), cudaMemcpyHostToDevice))

  /*****  Gauss-Seidel iterations *****/
  int iter = 0;            /* Current iteration number */
  double error = INFINITY; /* Current error */
  int has_not_converged = 1;

  // FREEZING CONTACTS
  if (iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
    printf("freezing contacts not supported yet\n");

    /***********************/
    /* NO FREEZING CONTACT */
    /***********************/
  } else {
    double light_error_sum = 0.;
    double* d_sumP2 = NULL;
    double* d_sumErr2 = NULL;

    double alpha = 1.0;
    double beta = 0.0;
    int n_rows;
    int start_line, end_line;
    int current_nnz;
    int cumul_nnz = 0;

    CHECK_CUDA(cudaMalloc(&d_sumP2, sizeof(double)))
    CHECK_CUDA(cudaMalloc(&d_sumErr2, sizeof(double)))

    CHECK_CUDA(cudaMemset(d_sumP2, 0, sizeof(double)))
    CHECK_CUDA(cudaMemset(d_sumErr2, 0, sizeof(double)))

    while ((iter < itermax) && has_not_converged) {
      // Set local problems q to q
      /* CHECK_CUDA(
          cudaMemcpy(d_q_new, d_q_const, (2 * nc) * sizeof(double), cudaMemcpyDeviceToDevice))
       */
      CHECK_CUDA(cudaMemset(d_sumP2, 0, sizeof(double)))
      CHECK_CUDA(cudaMemset(d_sumErr2, 0, sizeof(double)))

      cumul_nnz = 0;

      for (size_t color = 0; color < n_colors; color++) {
        int rangeSize = sum_sizes[color + 1] - sum_sizes[color];
        int threadsPerBlock = 256;
        int blocks = (rangeSize + threadsPerBlock - 1) / threadsPerBlock;
        size_t shmemSize = 2 * threadsPerBlock * sizeof(double);

        start_line = sum_sizes[color];
        end_line = sum_sizes[color + 1];
        n_rows = end_line - start_line;
        current_nnz =
            SBM_permuted->index1_data[end_line] - SBM_permuted->index1_data[start_line];

        CHECK_CUSPARSE(cusparseDbsrmv(
            handle, CUSPARSE_DIRECTION_COLUMN, CUSPARSE_OPERATION_NON_TRANSPOSE, n_rows,
            SBM_permuted->blocknumber1, current_nnz, &alpha, descr,
            &d_all_Values[4 * cumul_nnz], &d_all_RowOffsets[start_line + color],
            &d_all_ColIndices[cumul_nnz], 2, d_z, &beta, &d_q_new[2 * start_line]))

        fc2d_nsgs_local_solve_kernel_range_reduce_2<<<blocks, threadsPerBlock, shmemSize>>>(
            d_diagonal_blocks, d_determinants, d_q_new, d_q_const, d_mu, d_z, d_sumP2,
            d_sumErr2, sum_sizes[color], sum_sizes[color + 1]);

        cumul_nnz += current_nnz;
      }

      CHECK_CUDA(
          cudaMemcpy(&light_error_sum, d_sumErr2, sizeof(double), cudaMemcpyDeviceToHost))
      CHECK_CUDA(cudaMemcpy(&norm_r, d_sumP2, sizeof(double), cudaMemcpyDeviceToHost))

      error = sqrt(light_error_sum);
      norm_r = sqrt(norm_r);
      if (fabs(norm_r) > DBL_EPSILON) error /= norm_r;

      if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
          SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT) {
        has_not_converged = determine_convergence(error, tolerance, iter, options);
      } else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                 SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
        has_not_converged = determine_convergence_with_full_final(
            problem, options, z, d_z, w, &tolerance, norm_q, error, iter);
        problem->M->storageType = NM_SPARSE_BLOCK;
      }

      ++iter;
    }

    CHECK_CUDA(cudaFree(d_sumP2))
    CHECK_CUDA(cudaFree(d_sumErr2))
  }

  // Copy z back to host
  CHECK_CUDA(cudaMemcpy(z, d_z, (2 * nc) * sizeof(double), cudaMemcpyDeviceToHost))

  /* Full criterium */
  if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
      SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
    error = calculateFullErrorFinal(problem, options, z, w, tolerance, norm_q);
    problem->M->storageType = NM_SPARSE_BLOCK;

    has_not_converged =
        determine_convergence(error, dparam[SICONOS_DPARAM_TOL], iter, options);
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

  free(sum_sizes);
  free(inv_permutation);

  CHECK_CUSPARSE(cusparseDestroy(handle))

  /* Restore problem before permutation */
  problem->M->matrix1 = old_matrix1;
  problem->q = old_q;
  problem->mu = old_mu;

  if (SBM_problem != NULL) {
    SBM_free(SBM_problem, SBM_FREE_ALL);
    problem->M->matrix1 = NULL;
  }
  SBM_col_permuted =
      SBM_free(SBM_col_permuted, SBM_FREE_KEEP_BLOCKS);  // do not free blocks on this one
  SBM_permuted = SBM_free(SBM_permuted, SBM_FREE_ALL);  // free blocks because they were copied

  free(q_permuted);
  free(mu_permuted);

  // Device memory deallocation
  CHECK_CUDA(cudaFree(d_all_RowOffsets))
  CHECK_CUDA(cudaFree(d_all_ColIndices))
  CHECK_CUDA(cudaFree(d_all_Values))
  CHECK_CUDA(cudaFree(d_q_new))
  CHECK_CUDA(cudaFree(d_z))
  CHECK_CUDA(cudaFree(d_diagonal_blocks))
  CHECK_CUDA(cudaFree(d_determinants))
  CHECK_CUDA(cudaFree(d_mu))
  CHECK_CUDA(cudaFree(d_q_const))

  problem->M->storageType = old_storageType;
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers FC2D_NSGS_GRAPH_PERMUT_CUDA_BLOCKLEGACY in the global solver registry,
 * enabling:
 * - Dynamic solver lookup by ID
 * - Runtime solver introspection
 * - Elimination of giant switch statements in drivers
 */

static int fc2d_nsgs_graph_permut_cuda_blocklegacy_init_wrap(void* problem,
                                                             SolverOptions* options) {
  fc2d_nsgs_set_default(options);
  return NUMERICS_OK;
}

static int fc2d_nsgs_graph_permut_cuda_blocklegacy_solve_wrap(void* problem, double* reaction,
                                                              double* velocity,
                                                              SolverOptions* options) {
  int info = NUMERICS_OK;
  fc2d_nsgs_graph_permut_cuda_blocklegacy((FrictionContactProblem*)problem, reaction, velocity,
                                          &info, options);
  return info;
}

static void fc2d_nsgs_graph_permut_cuda_blocklegacy_free_wrap(void* problem,
                                                              SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_FRICTION_2D_NSGS_GRAPH_PERMUT_CUDA_BLOCKLEGACY,
                "SICONOS_FRICTION_2D_NSGS_GRAPH_PERMUT_CUDA_BLOCKLEGACY",
                "GPU implementation of FC2D_NSGS (uses legacy API for block format)",
                fc2d_nsgs_graph_permut_cuda_blocklegacy_init_wrap,
                fc2d_nsgs_graph_permut_cuda_blocklegacy_solve_wrap,
                fc2d_nsgs_graph_permut_cuda_blocklegacy_free_wrap, NULL, /* error function */
                fc2d_nsgs_set_default,                                   /* set_default */
                1000,                                                    /* default_max_iter */
                1e-4,                                                    /* default_tol */
                0);                                                      /* is_local_solver */