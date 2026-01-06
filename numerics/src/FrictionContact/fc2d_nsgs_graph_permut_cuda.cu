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
#include "graph_tools.h"
#include "numerics_verbose.h"  // for numerics_printf, verbose
#include "op3x3.h"

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
      return;                                                              \
    }                                                                      \
  }

#define CHECK_CUSPARSE(func)                                                   \
  {                                                                            \
    cusparseStatus_t status = (func);                                          \
    if (status != CUSPARSE_STATUS_SUCCESS) {                                   \
      printf("CUSPARSE API failed at line %d with error: %s (%d)\n", __LINE__, \
             cusparseGetErrorString(status), status);                          \
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

/* Set diagonal blocks to 0.
 * We could also try to completely remove them.
 */
static void remove_diagonal_blocks(FrictionContactProblem* problem) {
  unsigned int nc = problem->numberOfContacts;
  double* block;
  int diagPos;

  for (unsigned int i = 0; i < nc; ++i) {
    diagPos = SBM_diagonal_block_index(problem->M->matrix1, i);
    block = problem->M->matrix1->block[diagPos];
    for (int j = 0; j < 4; j++) block[j] = 0.0;
  }
}

/* CUDA kernel to solve local problems,
 * compute and aggregate the errors.
 */
__global__ void fc2d_nsgs_local_solve_kernel_range_reduce(
    const double* W,   // 4*n
    const double* D,   // size: n
    const double* q,   // 2*n
    const double* mu,  // size: n
    double* P,         // 2*n (updated in place)
    double* sumP2,     // scalar (device pointer)
    double* sumErr2,   // scalar (device pointer)
    const int k1, const int k2) {
  extern __shared__ double shmem[];
  double* shP2 = shmem;                 // blockDim.x
  double* shErr2 = shmem + blockDim.x;  // blockDim.x

  int tid = threadIdx.x;
  int local = blockIdx.x * blockDim.x + tid;
  int i = k1 + local;

  double localP2 = 0.0;
  double localErr2 = 0.0;

  if (i < k2) {
    const double Di = D[i];
    const double mui = mu[i];
    const double* Wi = &W[4 * i];
    const double* qi = &q[2 * i];
    double* Pi = &P[2 * i];

    // Save old P
    double Pold0 = Pi[0];
    double Pold1 = Pi[1];

    const double Wnn = Wi[0];
    const double Wtn = Wi[1];
    const double Wnt = Wi[2];
    const double Wtt = Wi[3];

    if (qi[0] > 0.0) {
      Pi[0] = 0.0;
      Pi[1] = 0.0;
    } else {
      Pi[0] = -(Wtt * qi[0] - Wnt * qi[1]) / Di;
      Pi[1] = -(-Wtn * qi[0] + Wnn * qi[1]) / Di;

      double muPn = mui * Pi[0];

      if (fabs(Pi[1]) > muPn) {
        if (Pi[1] + muPn < 0.0) {
          Pi[0] = -qi[0] / (Wnn - mui * Wnt);
          Pi[1] = -mui * Pi[0];
        } else {
          Pi[0] = -qi[0] / (Wnn + mui * Wnt);
          Pi[1] = mui * Pi[0];
        }
      }
    }

    // Accumulate local contributions
    double dP0 = Pi[0] - Pold0;
    double dP1 = Pi[1] - Pold1;

    localP2 = Pi[0] * Pi[0] + Pi[1] * Pi[1];
    localErr2 = dP0 * dP0 + dP1 * dP1;
  }

  // Store into shared memory
  shP2[tid] = localP2;
  shErr2[tid] = localErr2;
  __syncthreads();

  // Block reduction
  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (tid < s) {
      shP2[tid] += shP2[tid + s];
      shErr2[tid] += shErr2[tid + s];
    }
    __syncthreads();
  }

  // One atomic per block
  if (tid == 0) {
    atomicAdd(sumP2, shP2[0]);
    atomicAdd(sumErr2, shErr2[0]);
  }
}

void fc2d_nsgs_graph_permut_cuda(FrictionContactProblem* problem, double* z, double* w,
                                 int* info, SolverOptions* options) {
  /* Notes:
     - we suppose that the trivial solution case has been checked before,
     and that all inputs differs from NULL since this function is
     supposed to be called from lcp_driver_global().
  */

  // double time = omp_get_wtime();

  // Get solver parameters
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  unsigned int nc = problem->numberOfContacts;
  double norm_q = cblas_dnrm2(nc * 2, problem->q, 1);
  double norm_r = 0.0;
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = options->dparam[SICONOS_DPARAM_TOL];

  // time = omp_get_wtime() - time;
  // printf("time getting parameters: %es\n", time);
  // time = omp_get_wtime();

  // Initialize cusparse
  cusparseHandle_t handle = NULL;
  CHECK_CUSPARSE(cusparseCreate(&handle))
  // int cusparse_version = -1;
  // CHECK_CUSPARSE( cusparseGetVersion(handle, &cusparse_version) )
  // printf("cusparse version: %d\n", cusparse_version);

  // time = omp_get_wtime() - time;
  // printf("time initializing cusparse handle: %es\n", time);
  // time = omp_get_wtime();

  // Start coloring and permutation
  size_t n_colors = 0;
  size_t* sum_sizes = NULL;
  size_t* inv_permutation = (size_t*)malloc(nc * sizeof(size_t));

  color_graph_block_permut(nc, problem->M, &n_colors, &sum_sizes, inv_permutation);

  // time = omp_get_wtime() - time;
  // printf("time coloring: %es\n", time);

  // printf("Number of colors = %d\n", n_colors);

  // time = omp_get_wtime();

  SparseBlockStructuredMatrix* SBM_col_permuted = SBM_new();
  SparseBlockStructuredMatrix* SBM_permuted = SBM_new();
  unsigned int* rowIndex = (unsigned int*)malloc(nc * sizeof(unsigned int));

  for (unsigned int i = 0; i < nc; i++) rowIndex[inv_permutation[i]] = i;

  SBM_column_permutation(rowIndex, problem->M->matrix1, SBM_col_permuted);
  SBM_row_permutation_copy(inv_permutation, SBM_col_permuted, SBM_permuted);
  free(rowIndex);

  SparseBlockStructuredMatrix* old_matrix1 = problem->M->matrix1;
  problem->M->matrix1 = SBM_permuted;

  // time = omp_get_wtime() - time;
  // printf("time permutating: %e\n", time);

  // time = omp_get_wtime();

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

  // time = omp_get_wtime() - time;
  // printf("time extracting diagonal + determinants: %e\n", time);
  // time = omp_get_wtime();

  // This only sets elements to 0, but I think they will still be present (as 0s) in the CSR
  // conversion, unless the conversion function checks for 0 elements
  remove_diagonal_blocks(problem);  // set diagonal blocks of problem->M->matrix1 to 0

  CSparseMatrix* SBM_permuted_csr_nodiag = (CSparseMatrix*)malloc(sizeof(CSparseMatrix));
  size_t res = SBM_to_sparse_init_memory(problem->M->matrix1, SBM_permuted_csr_nodiag);
  res = SBM_to_sparse(problem->M->matrix1, SBM_permuted_csr_nodiag);
  /* if (res != 0) {
    printf("Error converting SBM to CSparse.\n");
    return;
  } else {
    printf("Conversion SBM->CSparse done\n");
  } */

  /*
    // Maximum number of elements in a line
    size_t max_nnz_row = 0;
    size_t nnz_row;
    for (size_t row = 0; row < 2 * nc; row++) {
      nnz_row = SBM_permuted_csr_nodiag->p[row + 1] - SBM_permuted_csr_nodiag->p[row];
      if (nnz_row > max_nnz_row) max_nnz_row = nnz_row;
    }
    printf("maximum number of elements in a row = %d\n", max_nnz_row);

    // Create sliced ELLpack matrices
    int64_t sliceSize = 8; // number of lines in a slice
    int64_t n_full_slices = (2 * nc) / sliceSize;
    int64_t last_slice_size = (2 * nc) % sliceSize;
    int64_t n_slices = n_full_slices + ((last_slice_size == 0) ? 0 : 1);

    int64_t nnz_sliced = 0;
    int64_t sellValuesSize = 0;


    for (unsigned int slice = 0; slice < n_full_slices; slice++) {
      // Get maximum number of elements in the slice
      size_t max = 0;
      size_t n_el;

      // Find max number of nonzero elements in a row for this slice
      for (unsigned int row = sliceSize * slice; row < sliceSize * (slice + 1); row++) {
        nnz_row = SBM_permuted_csr_nodiag->p[row + 1] - SBM_permuted_csr_nodiag->p[row];
        nnz_sliced += nnz_row;
        if (max < nnz_row) max = nnz_row;
      }

      sellValuesSize += max * sliceSize
    }

    // Last slice, possibly
    size_t max = 0;
    size_t n_el;

    for (unsigned int row = sliceSize * slice; row < sliceSize * (slice + 1); row++) {
        nnz_row = SBM_permuted_csr_nodiag->p[row + 1] - SBM_permuted_csr_nodiag->p[row];
        nnz_sliced += nnz_row;
        if (max < nnz_row) max = nnz_row;
      }

      sellValuesSize += max * sliceSize */

  // Create an array storing, for each color, the row offsets of the corresponding CSR
  // submatrix
  int64_t* h_all_RowOffsets = (int64_t*)malloc((2 * nc + n_colors) * sizeof(int64_t));
  size_t k = 0;
  for (unsigned int color = 0; color < n_colors; color++) {
    size_t start_line = 2 * sum_sizes[color];
    size_t end_line = 2 * sum_sizes[color + 1];

    for (unsigned int row = start_line; row < end_line + 1; row++) {
      h_all_RowOffsets[k] =
          SBM_permuted_csr_nodiag->p[row] - SBM_permuted_csr_nodiag->p[start_line];
      k++;
    }
  }

  // Copy this array to device
  int64_t* d_all_RowOffsets = NULL;
  CHECK_CUDA(cudaMalloc(&d_all_RowOffsets, (2 * nc + n_colors) * sizeof(int64_t)))
  CHECK_CUDA(cudaMemcpy(d_all_RowOffsets, h_all_RowOffsets,
                        (2 * nc + n_colors) * sizeof(int64_t), cudaMemcpyHostToDevice))
  free(h_all_RowOffsets);

  int64_t* d_all_ColIndices = NULL;
  double* d_all_Values = NULL;
  size_t nnz = SBM_permuted_csr_nodiag->p[2 * nc];

  CHECK_CUDA(cudaMalloc(&d_all_ColIndices, nnz * sizeof(int64_t)))
  CHECK_CUDA(cudaMalloc(&d_all_Values, nnz * sizeof(double)))

  CHECK_CUDA(cudaMemcpy(d_all_ColIndices, SBM_permuted_csr_nodiag->i, nnz * sizeof(int64_t),
                        cudaMemcpyHostToDevice))
  CHECK_CUDA(cudaMemcpy(d_all_Values, SBM_permuted_csr_nodiag->x, nnz * sizeof(double),
                        cudaMemcpyHostToDevice))

  cusparseConstSpMatDescr_t* all_cusparse_csrmat =
      (cusparseConstSpMatDescr_t*)malloc(n_colors * sizeof(cusparseConstSpMatDescr_t));

  double alpha = 1.0;
  double beta = 1.0;
  size_t n_rows;
  size_t start_line, end_line;
  size_t current_nnz;
  size_t cumul_nnz = 0;
  size_t cumul_n_rows = 0;
  for (unsigned int color = 0; color < n_colors; color++) {
    start_line = 2 * sum_sizes[color];
    end_line = 2 * sum_sizes[color + 1];
    n_rows = end_line - start_line;
    cumul_n_rows += n_rows;
    current_nnz =
        SBM_permuted_csr_nodiag->p[end_line] - SBM_permuted_csr_nodiag->p[start_line];

    // printf("%d\n", h_all_RowOffsets[start_line + color]);

    CHECK_CUSPARSE(cusparseCreateConstCsr(
        &all_cusparse_csrmat[color], n_rows, 2 * nc, current_nnz,
        &d_all_RowOffsets[start_line + color], &d_all_ColIndices[cumul_nnz],
        &d_all_Values[cumul_nnz], CUSPARSE_INDEX_64I, CUSPARSE_INDEX_64I,
        CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F))

    cumul_nnz += current_nnz;
  }

  // Create vector for z, the reaction, which is the input of the product
  double* d_z = NULL;
  CHECK_CUDA(cudaMalloc(&d_z, (2 * nc) * sizeof(double)))
  CHECK_CUDA(cudaMemcpy(d_z, z, (2 * nc) * sizeof(double), cudaMemcpyHostToDevice))
  cusparseDnVecDescr_t vec_z;
  CHECK_CUSPARSE(cusparseCreateDnVec(&vec_z, 2 * nc, d_z, CUDA_R_64F))
  // z vector created

  // Create the res vectors, which are the q's of the local problems
  double* d_q_new = NULL;
  CHECK_CUDA(cudaMalloc(&d_q_new, (2 * nc) * sizeof(double)))
  CHECK_CUDA(
      cudaMemcpy(d_q_new, problem->q, (2 * nc) * sizeof(double), cudaMemcpyHostToDevice))

  cusparseDnVecDescr_t* all_q_new =
      (cusparseDnVecDescr_t*)malloc(n_colors * sizeof(cusparseDnVecDescr_t));

  for (unsigned int color = 0; color < n_colors; color++) {
    CHECK_CUSPARSE(cusparseCreateDnVec(&all_q_new[color],
                                       2 * (sum_sizes[color + 1] - sum_sizes[color]),
                                       &d_q_new[2 * sum_sizes[color]], CUDA_R_64F))
  }

  void** d_all_Buffer = (void**)malloc(n_colors * sizeof(void*));
  size_t* all_Buffer_Size = (size_t*)malloc(n_colors * sizeof(size_t));

  alpha = 1.0;
  beta = 1.0;
  for (unsigned int color = 0; color < n_colors; color++) {
    CHECK_CUSPARSE(cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha,
                                           all_cusparse_csrmat[color], vec_z, &beta,
                                           all_q_new[color], CUDA_R_64F,
                                           CUSPARSE_SPMV_ALG_DEFAULT, &all_Buffer_Size[color]))
    CHECK_CUDA(cudaMalloc(&d_all_Buffer[color], all_Buffer_Size[color]))

    CHECK_CUSPARSE(cusparseSpMV_preprocess(
        handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, all_cusparse_csrmat[color], vec_z,
        &beta, all_q_new[color], CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, d_all_Buffer[color]))
  }

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

  // time = omp_get_wtime() - time;
  // printf("time initializing cusparse stuff: %es\n", time);

  /*****  Gauss-Seidel iterations *****/
  int iter = 0;            /* Current iteration number */
  double error = INFINITY; /* Current error */
  int has_not_converged = 1;

  // double time = omp_get_wtime();

  // FREEZING CONTACTS
  if (iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
    printf("freezing contacts not supported yet\n");
    return;
    /* unsigned int contact;
    unsigned int pos;
    double light_error_sum;
    double light_error_2;
    double localreaction[2];
    unsigned int number_of_freezed_contact = 0;
    double tmp_criteria1, tmp_criteria2;

    freeze_contacts = f2d_nsgs_allocate_freezing_contacts(problem, options);

#pragma omp parallel default(none) private(pos, local_problem, localreaction, light_error_2, \
                                               index1_data, index2_data)                     \
    shared(problem, diagonal_blocks, diagonal_block_determinant, z, sum_sizes, iter)         \
    shared(iparam, light_error_sum, n_colors, norm_r, nc, error, options, tolerance,         \
               has_not_converged, norm_q, w, itermax)                                        \
    shared(tmp_criteria1, tmp_criteria2, freeze_contacts, number_of_freezed_contact,         \
               blocks_contiguous)
    {
      local_problem =
          (LinearComplementarityProblem*)malloc(sizeof(LinearComplementarityProblem));
      local_problem->M = NM_new();
      local_problem->M->storageType = NM_DENSE;
      local_problem->M->size0 = 2;
      local_problem->M->size1 = 2;
      local_problem->q = (double*)malloc(2 * sizeof(double));

      while ((iter < itermax) && has_not_converged) {
        // light_error_sum = 0.0;
        light_error_2 = 0.0;

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
          for (unsigned int permuted_contact = sum_sizes[color];
               permuted_contact < sum_sizes[color + 1]; permuted_contact++) {
            if (freeze_contacts[permuted_contact] > 0) {
              freeze_contacts[permuted_contact] -= 1;
              continue;
            }

            pos = 2 * permuted_contact;
            localreaction[0] = z[pos];
            localreaction[1] = z[pos + 1];

            fc2d_nsgs_buildLocalProblem_parallel(permuted_contact, problem, blocks_contiguous,
                                                 diagonal_blocks, index1_data, index2_data,
                                                 local_problem, z);

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
              freeze_contacts[permuted_contact] =
                  iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT];
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

          if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
              SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT) {
            error = calculateLightError(light_error_sum, nc, z, &norm_r);
            has_not_converged = determine_convergence(error, tolerance, iter, options);
          } else if (iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] ==
                     SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL) {
            error = calculateLightError(light_error_sum, nc, z, &norm_r);
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
    }  // end parallel region */

    /***********************/
    /* NO FREEZING CONTACT */
    /***********************/
  } else {
    double light_error_sum = 0.;
    double* d_sumP2 = NULL;
    double* d_sumErr2 = NULL;

    CHECK_CUDA(cudaMalloc(&d_sumP2, sizeof(double)))
    CHECK_CUDA(cudaMalloc(&d_sumErr2, sizeof(double)))

    CHECK_CUDA(cudaMemset(d_sumP2, 0, sizeof(double)))
    CHECK_CUDA(cudaMemset(d_sumErr2, 0, sizeof(double)))

    while ((iter < itermax) && has_not_converged) {
      // Set local problems q to q
      CHECK_CUDA(
          cudaMemcpy(d_q_new, d_q_const, (2 * nc) * sizeof(double), cudaMemcpyDeviceToDevice))
      CHECK_CUDA(cudaMemset(d_sumP2, 0, sizeof(double)))
      CHECK_CUDA(cudaMemset(d_sumErr2, 0, sizeof(double)))

      for (size_t color = 0; color < n_colors; color++) {
        int rangeSize = sum_sizes[color + 1] - sum_sizes[color];
        int threadsPerBlock = 256;
        int blocks = (rangeSize + threadsPerBlock - 1) / threadsPerBlock;
        size_t shmemSize = 2 * threadsPerBlock * sizeof(double);

        // Build local problems in parallel
        CHECK_CUSPARSE(cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha,
                                    all_cusparse_csrmat[color], vec_z, &beta, all_q_new[color],
                                    CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT,
                                    d_all_Buffer[color]))

        fc2d_nsgs_local_solve_kernel_range_reduce<<<blocks, threadsPerBlock, shmemSize>>>(
            d_diagonal_blocks, d_determinants, d_q_new, d_mu, d_z, d_sumP2, d_sumErr2,
            sum_sizes[color], sum_sizes[color + 1]);
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
            problem, options, z, w, &tolerance, norm_q, error, iter);
      }

      ++iter;
    }

    CHECK_CUDA(cudaFree(d_sumP2))
    CHECK_CUDA(cudaFree(d_sumErr2))
  }

  /* time = omp_get_wtime() - time;
  printf("time in loop: %es\n", time); */

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

  // Destroy matrix/vector descriptors
  for (unsigned int color = 0; color < n_colors; color++) {
    CHECK_CUSPARSE(cusparseDestroySpMat(all_cusparse_csrmat[color]))
    CHECK_CUSPARSE(cusparseDestroyDnVec(all_q_new[color]))
    CHECK_CUDA(cudaFree(d_all_Buffer[color]))
  }
  free(all_cusparse_csrmat);
  free(all_q_new);
  free(d_all_Buffer);
  free(all_Buffer_Size);

  CHECK_CUSPARSE(cusparseDestroyDnVec(vec_z))
  CHECK_CUSPARSE(cusparseDestroy(handle))

  /* Restore problem before permutation */
  problem->M->matrix1 = old_matrix1;
  problem->q = old_q;
  problem->mu = old_mu;

  free(SBM_col_permuted);
  free(SBM_permuted);
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
}
