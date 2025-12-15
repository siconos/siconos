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

/* Only support SBM format for now */
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

static void shuffle(unsigned int size, unsigned int* randnum)  // size is the given range
{
  unsigned int swap, randindex;
  for (unsigned i = 0; i < size; ++i) {
    swap = randnum[i];
    randindex = rand() % size;
    randnum[i] = randnum[randindex];
    randnum[randindex] = swap;
  }
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
  *norm_r = cblas_dnrm2(nc * 2, reaction, 1);
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
  unsigned int nc = problem->numberOfContacts;
  if (options->iparam[SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] > 0) {
    fcontacts = (unsigned int*)malloc(nc * sizeof(unsigned int));
    for (unsigned int i = 0; i < nc; ++i) {
      fcontacts[i] = 0;
    }
  }
  return fcontacts;
}

static void remove_diagonal_blocks(FrictionContactProblem *problem)
{
    unsigned int nc = problem->numberOfContacts;
    double *block;
    int diagPos;

    for (unsigned int i = 0; i < nc; ++i)
    {
        diagPos = SBM_diagonal_block_index(problem->M->matrix1, i);
        block = problem->M->matrix1->block[diagPos];
        for (int j = 0; j < 4; j++)
            block[j] = 0.0;
    }
}

// CUDA kernel for vector addition
__global__ void vecAdd(float* A, float* B, float* C, int N) {
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  if (i < N) {
    C[i] = A[i] + B[i];
  }
}

// Device function to solve local problem
__device__ __forceinline__ void fc2d_nsgs_local_solve_device(
    const double* W,
    double D,
    const double* q,
    double mu,
    double* P)
{
    // | Wnn Wnt |
    // | Wtn Wtt |

    const double Wnn = W[0];
    const double Wtn = W[1];
    const double Wnt = W[2];
    const double Wtt = W[3];

    if (q[0] > 0.0) {
        P[0] = 0.0;
        P[1] = 0.0;
    } else {
        // solve W P + q = 0
        P[0] = -(Wtt * q[0] - Wnt * q[1]) / D;
        P[1] = -(-Wtn * q[0] + Wnn * q[1]) / D;

        const double muPn = mu * P[0];

        // outside cone
        if (fabs(P[1]) > muPn) {
            if (P[1] + muPn < 0.0) {
                P[0] = -q[0] / (Wnn - mu * Wnt);
                P[1] = -mu * P[0];
            } else {
                P[0] = -q[0] / (Wnn + mu * Wnt);
                P[1] =  mu * P[0];
            }
        }
    }
}

// CUDA kernel to solve local problem
__global__ void fc2d_nsgs_local_solve_kernel(
    const double* W,   // size: 4 * n
    const double* D,   // size: n
    const double* q,   // size: 2 * n
    const double* mu,  // size: n
    double* P,         // size: 2 * n
    int k1,
    int k2)
{
    int local = blockIdx.x * blockDim.x + threadIdx.x;
    int i = k1 + local;
    if (i >= k2) return;

    fc2d_nsgs_local_solve_device(
        &W[4 * i],
        D[i],
        &q[2 * i],
        mu[i],
        &P[2 * i]
    );
}


__global__
void fc2d_nsgs_local_solve_kernel_range_reduce(
    const double* W,      // 4*n
    const double* D,      // size: n
    const double* q,      // 2*n
    const double* mu,     // size: n
    double* P,            // 2*n (updated in place)
    double* sumP2,        // scalar (device pointer)
    double* sumErr2,      // scalar (device pointer)
    int k1,
    int k2)
{
    extern __shared__ double shmem[];
    double* shP2   = shmem;                 // blockDim.x
    double* shErr2 = shmem + blockDim.x;    // blockDim.x

    int tid   = threadIdx.x;
    int local = blockIdx.x * blockDim.x + tid;
    int i     = k1 + local;

    double localP2   = 0.0;
    double localErr2 = 0.0;

    if (i < k2) {
        const double Di = D[i];
        const double mui = mu[i];
        const double* Wi = &W[4 * i];
        const double* qi = &q[2 * i];
        double* Pi       = &P[2 * i];

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
                    Pi[1] =  mui * Pi[0];
                }
            }
        }

        // Accumulate local contributions
        double dP0 = Pi[0] - Pold0;
        double dP1 = Pi[1] - Pold1;

        localP2   = Pi[0] * Pi[0] + Pi[1] * Pi[1];
        localErr2 = dP0 * dP0 + dP1 * dP1;
    }

    // Store into shared memory
    shP2[tid]   = localP2;
    shErr2[tid] = localErr2;
    __syncthreads();

    // Block reduction
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            shP2[tid]   += shP2[tid + s];
            shErr2[tid] += shErr2[tid + s];
        }
        __syncthreads();
    }

    // One atomic per block
    if (tid == 0) {
        atomicAdd(sumP2,   shP2[0]);
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

  // Get solver parameters
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  unsigned int nc = problem->numberOfContacts;
  double norm_q = cblas_dnrm2(nc * 2, problem->q, 1);
  double norm_r = 0.0;
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tolerance = options->dparam[SICONOS_DPARAM_TOL];

  // Initialize cusparse
  cusparseHandle_t handle;
  CHECK_CUSPARSE( cusparseCreate(&handle) );
  int cusparse_version = -1;
  CHECK_CUSPARSE( cusparseGetVersion(handle, &cusparse_version) )
  printf("cusparse version: %d\n", cusparse_version);

  /* *********** */
  /* vecAdd code */
  /* *********** */

  /*
  int N = 1024;
  size_t size = N * sizeof(float);

  // Host input vectors
  float *h_A, *h_B;
  // Host output vector
  float* h_C;
  // Device input vectors
  float *d_A, *d_B;
  // Device output vector
  float* d_C;

  // Allocate memory for each vector on host
  h_A = (float*)malloc(size);
  h_B = (float*)malloc(size);
  h_C = (float*)malloc(size);

  // Allocate memory for each vector on GPU
  cudaMalloc(&d_A, size);
  cudaMalloc(&d_B, size);
  cudaMalloc(&d_C, size);

  // Initialize vectors on host
  for (int i = 0; i < N; i++) {
    h_A[i] = sin(i) * sin(i);
    h_B[i] = cos(i) * cos(i);
  }

  // Copy host vectors to device
  cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice);

  // Execute the kernel
  int threadsPerBlock = 256;
  int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
  vecAdd<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, N);

  // Copy array back to host
  cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);

  // Sum up vector c and print result divided by N, this should equal 1 within error
  float sum = 0;
  for (int i = 0; i < N; i++) {
    sum += h_C[i];
  }

  printf("Final result: %e\n", sum / N);

  // Release device memory
  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_C);

  // Release host memory
  free(h_A);
  free(h_B);
  free(h_C); */

  /* ********************* */
  /* END OF VECADD EXAMPLE */
  /* ********************* */

  /* Example of matrix vector product to test stuff */
  /* double alpha = 1.0;
  double beta = 0.0;
  double* h_res = (double*)calloc(nc * 2, sizeof(double));
  NM_gemv(alpha, problem->M, problem->q, beta, h_res);
  for (int i = 0; i < 10; i++) {
    printf(" %e ", h_res[i]);
  }
  printf("\n"); */

  /* SpMV with coloring example */

  /* const int64_t A_num_rows      = 4;
  const int64_t A_num_cols      = 4;
  const int64_t A_nnz           = 9;
  int64_t       hA_csrOffsets[] = { 0, 3, 4, 7, 9 };
  int64_t       hA_columns[]    = { 0, 2, 3, 1, 0, 2, 3, 1, 3 };
  double     hA_values[]     = { 1.0, 2.0, 3.0, 4.0, 5.0,
                                6.0, 7.0, 8.0, 9.0 };
  double     hX[]            = { 1.0, 2.0, 3.0, 4.0 };
  double     hY[]            = { 0.0, 0.0, 0.0, 0.0 };
  double     hY_result[]     = { 19.0, 8.0, 51.0, 52.0 };
  int64_t       nb_of_colors    = 2;
  int64_t        partitions[]    = {0, 1, 2};
  //--------------------------------------------------------------------------
  // Device memory management
  cusparseHandle_t     handle = NULL;
  CHECK_CUSPARSE( cusparseCreate(&handle) )
  int64_t   *dA_csrOffsets, *dA_columns;
  double *dA_values, *dX, *dY;
  CHECK_CUDA( cudaMalloc((void**) &dA_csrOffsets,
                          (A_num_rows + 1) * sizeof(int64_t)) )
  CHECK_CUDA( cudaMalloc((void**) &dA_columns, A_nnz * sizeof(int64_t))        )
  CHECK_CUDA( cudaMalloc((void**) &dA_values,  A_nnz * sizeof(double))      )
  CHECK_CUDA( cudaMalloc((void**) &dX,         A_num_cols * sizeof(double)) )
  CHECK_CUDA( cudaMalloc((void**) &dY,         A_num_rows * sizeof(double)) )

  CHECK_CUDA( cudaMemcpy(dA_csrOffsets, hA_csrOffsets,
                          (A_num_rows + 1) * sizeof(int64_t),
                          cudaMemcpyHostToDevice) )
  CHECK_CUDA( cudaMemcpy(dA_columns, hA_columns, A_nnz * sizeof(int64_t),
                          cudaMemcpyHostToDevice) )
  CHECK_CUDA( cudaMemcpy(dA_values, hA_values, A_nnz * sizeof(double),
                          cudaMemcpyHostToDevice) )
  CHECK_CUDA( cudaMemcpy(dX, hX, A_num_cols * sizeof(double),
                          cudaMemcpyHostToDevice) )
  CHECK_CUDA( cudaMemcpy(dY, hY, A_num_rows * sizeof(double),
                          cudaMemcpyHostToDevice) )

  // Create dense vector X
  cusparseDnVecDescr_t vecX;
  CHECK_CUSPARSE( cusparseCreateDnVec(&vecX, A_num_cols, dX, CUDA_R_64F) )

  int64_t *h_all_RowOffsets = (int64_t *)malloc((A_num_rows + nb_of_colors) * sizeof(int64_t));
  unsigned int k = 0;
  size_t start_line;
  size_t end_line;
  for (unsigned int color = 0; color < nb_of_colors; color++) {
    start_line = 2 * partitions[color];
    end_line = 2 * partitions[color + 1];

    for (unsigned int row = start_line; row < end_line + 1; row++) {
      h_all_RowOffsets[k] = hA_csrOffsets[row] - hA_csrOffsets[start_line];
      k++;
    }
  }

  for (int i = 0; i < A_num_rows + nb_of_colors; i++) printf(" %d ", h_all_RowOffsets[i]);
  printf("\n");

  int64_t *d_all_RowOffsets = NULL;
  CHECK_CUDA( cudaMalloc(&d_all_RowOffsets, (A_num_rows + nb_of_colors) * sizeof(int64_t)) )
  CHECK_CUDA( cudaMemcpy(d_all_RowOffsets, h_all_RowOffsets, (A_num_rows + nb_of_colors) * sizeof(int64_t),
             cudaMemcpyHostToDevice) )

  int64_t* d_all_ColIndices = NULL;
  double* d_all_Values = NULL;

  CHECK_CUDA( cudaMalloc(&d_all_ColIndices, A_nnz * sizeof(int64_t)) )
  CHECK_CUDA( cudaMalloc(&d_all_Values, A_nnz * sizeof(double)) )

  CHECK_CUDA( cudaMemcpy(d_all_ColIndices, hA_columns, A_nnz * sizeof(int64_t),
             cudaMemcpyHostToDevice) )
  CHECK_CUDA( cudaMemcpy(d_all_Values, hA_values, A_nnz * sizeof(double),
             cudaMemcpyHostToDevice) )

  cusparseConstSpMatDescr_t *all_cusparse_csrmat = (cusparseConstSpMatDescr_t *)malloc(nb_of_colors * sizeof(cusparseConstSpMatDescr_t));
  
  size_t n_rows;
  size_t current_nnz;
  size_t cumul_nnz = 0;
  size_t cumul_n_rows = 0;
  for (unsigned int color = 0; color < nb_of_colors; color++) {
    start_line = 2 * partitions[color];
    end_line = 2 * partitions[color + 1];
    n_rows = end_line - start_line;
    cumul_n_rows += n_rows;
    current_nnz = hA_csrOffsets[end_line] - hA_csrOffsets[start_line];

    printf("%d\n", h_all_RowOffsets[start_line + color]);

    CHECK_CUSPARSE( cusparseCreateConstCsr(&all_cusparse_csrmat[color], n_rows, A_num_cols, current_nnz, 
                                           &d_all_RowOffsets[start_line + color], &d_all_ColIndices[cumul_nnz], &d_all_Values[cumul_nnz], 
                                          CUSPARSE_INDEX_64I, CUSPARSE_INDEX_64I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) )

    cumul_nnz += current_nnz;
  }

  printf("cumul_n_rows = %d, n_rows = %d\n", cumul_n_rows, A_num_rows);
  printf("cumul_nnz = %d, nnz = %d\n", cumul_nnz, A_nnz);

  printf("matrices created\n");
  fflush(stdout);

  // Create a res vector for each color
  double *h_all_res = (double *)calloc(A_num_rows, sizeof(double));
  double *d_all_res = NULL;
  CHECK_CUDA( cudaMalloc(&d_all_res, (A_num_rows) * sizeof(double)) )
  CHECK_CUDA( cudaMemcpy(d_all_res, h_all_res, (A_num_rows) * sizeof(double), cudaMemcpyHostToDevice) )

  cusparseDnVecDescr_t *all_vec_res = (cusparseDnVecDescr_t *)malloc(nb_of_colors * sizeof(cusparseDnVecDescr_t));
  for (unsigned int color = 0; color < nb_of_colors; color++) {
    CHECK_CUSPARSE( cusparseCreateDnVec(&all_vec_res[color], 2 * (partitions[color + 1] - partitions[color]), &d_all_res[2 * partitions[color]], CUDA_R_64F) )
  }

  printf("vectors created\n");
  fflush(stdout);
  
  void **d_all_Buffer = (void **)malloc(nb_of_colors * sizeof(void *));
  size_t *all_Buffer_Size = (size_t *)malloc(nb_of_colors * sizeof(size_t));

  for (unsigned int color = 0; color < nb_of_colors; color++) {
    CHECK_CUSPARSE( cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, all_cusparse_csrmat[color],
                                          vecX, &beta, all_vec_res[color], CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, &all_Buffer_Size[color]) )
    CHECK_CUDA( cudaMalloc(&d_all_Buffer[color], all_Buffer_Size[color]) )

    CHECK_CUSPARSE( cusparseSpMV_preprocess(
                                handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha, all_cusparse_csrmat[color], vecX, &beta, all_vec_res[color], CUDA_R_64F,
                                CUSPARSE_SPMV_ALG_DEFAULT, d_all_Buffer[color]) )

    CHECK_CUSPARSE( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, all_cusparse_csrmat[color],
                                 vecX, &beta, all_vec_res[color], CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, d_all_Buffer[color]) )
  }
  

  CHECK_CUDA( cudaMemcpy(h_all_res, d_all_res, A_num_rows * sizeof(double), cudaMemcpyDeviceToHost) )

  int correct = 1;
  for (int i = 0; i < A_num_rows; i++) {
    printf("%e %e\n", h_all_res[i], hY_result[i]);
      if (fabs(h_all_res[i] - hY_result[i]) / fabs(hY_result[i]) > 1e-15) { // direct floating point comparison is not
        correct = 0;             // reliable
        break;
      }
  }
  if (correct)
      printf("Dummy multicolor test PASSED\n");
  else
      printf("Dummy multicolor test FAILED\n"); */
  //--------------------------------------------------------------------------
  // CUSPARSE APIs
  
  /* cusparseSpMatDescr_t matA;
  cusparseDnVecDescr_t vecX, vecY;
  void*                dBuffer    = NULL;
  size_t               bufferSize = 0;
  CHECK_CUSPARSE( cusparseCreate(&handle) )
  // Create sparse matrix A in CSR format
  CHECK_CUSPARSE( cusparseCreateCsr(&matA, A_num_rows, A_num_cols, A_nnz,
                                    dA_csrOffsets, dA_columns, dA_values,
                                    CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                    CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) )
  // Create dense vector X
  CHECK_CUSPARSE( cusparseCreateDnVec(&vecX, A_num_cols, dX, CUDA_R_64F) )
  // Create dense vector y
  CHECK_CUSPARSE( cusparseCreateDnVec(&vecY, A_num_rows, dY, CUDA_R_64F) )
  // allocate an external buffer if needed
  CHECK_CUSPARSE( cusparseSpMV_bufferSize(
                                handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
                                CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize) )
  CHECK_CUDA( cudaMalloc(&dBuffer, bufferSize) )

  // execute preprocess (optional)
  CHECK_CUSPARSE( cusparseSpMV_preprocess(
                                handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
                                CUSPARSE_SPMV_ALG_DEFAULT, dBuffer) )

  // execute SpMV
  CHECK_CUSPARSE( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
                                CUSPARSE_SPMV_ALG_DEFAULT, dBuffer) )

  // destroy matrix/vector descriptors
  CHECK_CUSPARSE( cusparseDestroySpMat(matA) )
  CHECK_CUSPARSE( cusparseDestroyDnVec(vecX) )
  CHECK_CUSPARSE( cusparseDestroyDnVec(vecY) )

  //--------------------------------------------------------------------------
  // device result check
  CHECK_CUDA( cudaMemcpy(hY, dY, A_num_rows * sizeof(double),
                          cudaMemcpyDeviceToHost) )
  int correct = 1;
  for (int i = 0; i < A_num_rows; i++) {
      printf("%e %e\n", hY[i], hY_result[i]);
      if (hY[i] != hY_result[i]) { // direct floating point comparison is not
          correct = 0;             // reliable
          break;
      }
  }
  if (correct)
      printf("spmv_csr_example test PASSED\n");
  else
      printf("spmv_csr_example test FAILED: wrong result\n");
  //--------------------------------------------------------------------------
  // device memory deallocation
  CHECK_CUDA( cudaFree(dBuffer) )
  CHECK_CUDA( cudaFree(dA_csrOffsets) )
  CHECK_CUDA( cudaFree(dA_columns) )
  CHECK_CUDA( cudaFree(dA_values) )
  CHECK_CUDA( cudaFree(dX) )
  CHECK_CUDA( cudaFree(dY) ) */

  /* *********************** */
  /* END OF SMPV CSR EXAMPLE */
  /* *********************** */

  /* ******** */
  /* CUSPARSE */
  /* ******** */

  /* CSparseMatrix *csparse_mat = (CSparseMatrix *)malloc(sizeof(CSparseMatrix));
  int res = SBM_to_sparse_init_memory(problem->M->matrix1, csparse_mat);
  res = SBM_to_sparse(problem->M->matrix1, csparse_mat);
  if (res != 0) {
    printf("Error converting SBM to CSparse.\n");
    return;
  }
  else {
    printf("Conversion SBM->CSparse done\n");
  }

  int64_t* d_RowOffsets = NULL;
  int64_t* d_ColIndices = NULL;
  double* d_Values = NULL;
  size_t nnz = csparse_mat->p[2 * nc];
  printf("nnz = %d, nbblocks * 4 = %d\n", nnz, problem->M->matrix1->nbblocks * 4);

  CHECK_CUDA( cudaMalloc(&d_RowOffsets, (2 * nc + 1) * sizeof(int64_t)) )
  CHECK_CUDA( cudaMalloc(&d_ColIndices, nnz * sizeof(int64_t)) )
  CHECK_CUDA( cudaMalloc(&d_Values, nnz * sizeof(double)) )

  CHECK_CUDA( cudaMemcpy(d_RowOffsets, csparse_mat->p, (2 * nc + 1) * sizeof(int64_t),
             cudaMemcpyHostToDevice) )
  CHECK_CUDA( cudaMemcpy(d_ColIndices, csparse_mat->i, nnz * sizeof(int64_t),
             cudaMemcpyHostToDevice) )
  CHECK_CUDA( cudaMemcpy(d_Values, csparse_mat->x, nnz * sizeof(double),
             cudaMemcpyHostToDevice) )

  cusparseConstSpMatDescr_t cusparse_csrmat;
  CHECK_CUSPARSE( cusparseCreateConstCsr(&cusparse_csrmat, 2 * nc, 2 * nc, nnz, d_RowOffsets, d_ColIndices, d_Values, CUSPARSE_INDEX_64I, CUSPARSE_INDEX_64I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) )
  printf("cusparse CSR matrix created.\n");

  double* d_q = NULL;
  CHECK_CUDA( cudaMalloc(&d_q, (2 * nc) * sizeof(double)) )
  CHECK_CUDA( cudaMemcpy(d_q, problem->q, (2 * nc) * sizeof(double), cudaMemcpyHostToDevice) )
  cusparseDnVecDescr_t vecq;
  CHECK_CUSPARSE( cusparseCreateDnVec(&vecq, 2 * nc, d_q, CUDA_R_64F) )

  double *h_res_2 = (double *)calloc(2 * nc, sizeof(double));
  double *d_res = NULL;
  CHECK_CUDA( cudaMalloc(&d_res, (2 * nc) * sizeof(double)) )
  CHECK_CUDA( cudaMemcpy(d_res, h_res_2, (2 * nc) * sizeof(double), cudaMemcpyHostToDevice) )
  cusparseDnVecDescr_t vecres;
  CHECK_CUSPARSE( cusparseCreateDnVec(&vecres, 2 * nc, d_res, CUDA_R_64F) )

  void *dBuffer_2 = NULL;
  size_t bufferSize_2 = 0;

  CHECK_CUSPARSE( cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, cusparse_csrmat,
                                          vecq, &beta, vecres, CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize_2) )
  CHECK_CUDA( cudaMalloc(&dBuffer_2, bufferSize_2) )

  CHECK_CUSPARSE( cusparseSpMV_preprocess(
                                handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha, cusparse_csrmat, vecq, &beta, vecres, CUDA_R_64F,
                                CUSPARSE_SPMV_ALG_DEFAULT, dBuffer_2) )

  CHECK_CUSPARSE( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, cusparse_csrmat,
                               vecq, &beta, vecres, CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer_2) )

  // CHECK_CUSPARSE( cusparseDestroySpMat(cusparse_csrmat) )
  // CHECK_CUSPARSE( cusparseDestroyDnVec(vecq) )
  // CHECK_CUSPARSE( cusparseDestroyDnVec(vecres) )

  // CHECK_CUSPARSE( cusparseDnVecGetValues(vecq, d_res) );
  CHECK_CUDA( cudaMemcpy(h_res_2, d_res, (2 * nc) * sizeof(double), cudaMemcpyDeviceToHost) )

  for (int i = 0; i < 10; i++) printf(" %e ", h_res_2[i]);
  printf("\n");

  printf("DBL_EPSILON = %e\n", DBL_EPSILON);
  correct = 1;
  for (int i = 0; i < 2 * nc; i++) {
      if (fabs(h_res[i] - h_res_2[i]) / fabs(h_res[i]) > 1e-10) { // direct floating point comparison is not
          printf("%d %e %e %e\n", i, h_res[i], h_res_2[i], fabs(h_res[i] - h_res_2[i]) / fabs(h_res[i]));
          correct = 0;             // reliable
          break;
      }
  }
  if (correct)
      printf("Conversion test PASSED\n");
  else
      printf("Conversion test FAILED\n"); */
  
  // TEST WITH COLORING + PERMUTATION

  // Start coloring and permutation 
  size_t n_colors = 0;
  size_t* sum_sizes = NULL;
  size_t* inv_permutation = (size_t*)malloc(nc * sizeof(size_t));

  color_graph_block_permut(nc, problem->M, &n_colors, &sum_sizes, inv_permutation);

  numerics_printf("-- FC2D NSGS CUDA - Graph colored with %zu colors\n", n_colors);

  SparseBlockStructuredMatrix* SBM_col_permuted = SBM_new();
  SparseBlockStructuredMatrix* SBM_permuted = SBM_new();
  unsigned int* rowIndex = (unsigned int*)malloc(nc * sizeof(unsigned int));

  for (unsigned int i = 0; i < nc; i++) rowIndex[inv_permutation[i]] = i;

  SBM_column_permutation(rowIndex, problem->M->matrix1, SBM_col_permuted);
  SBM_row_permutation_copy(inv_permutation, SBM_col_permuted, SBM_permuted);
  free(rowIndex);

  SparseBlockStructuredMatrix* old_matrix1 = problem->M->matrix1;
  problem->M->matrix1 = SBM_permuted;

  /* Get diagonal blocks and determinants */
  double* diagonal_blocks = fc2d_extract_diagonal_blocks(problem);
  double* diagonal_block_determinant = fc2d_nsgs_compute_local_problem_determinant(nc, diagonal_blocks);

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

  /* double alpha = 1.0;
  double beta = 0.0;
  double* h_res_permuted = (double*)calloc(nc * 2, sizeof(double));
  NM_gemv(alpha, problem->M, problem->q, beta, h_res_permuted);
  for (int i = 0; i < 10; i++) {
    printf(" %e ", h_res_permuted[i]);
  }
  printf("\n");

  // Conversion to CSR
  CSparseMatrix *SBM_permuted_csr = (CSparseMatrix *)malloc(sizeof(CSparseMatrix));
  int res = SBM_to_sparse_init_memory(SBM_permuted, SBM_permuted_csr);
  res = SBM_to_sparse(SBM_permuted, SBM_permuted_csr);
  if (res != 0) {
    printf("Error converting SBM to CSparse.\n");
    return;
  }
  else {
    printf("Conversion SBM->CSparse done\n");
  }
  // Conversion done

  // Create an array storing, for each color, the row offsets of the corresponding CSR submatrix 
  int64_t *h_all_RowOffsets_ = (int64_t *)malloc((2 * nc + n_colors) * sizeof(int64_t));
  int k = 0;
  for (unsigned int color = 0; color < n_colors; color++) {
    size_t start_line = 2 * sum_sizes[color];
    size_t end_line = 2 * sum_sizes[color + 1];

    for (unsigned int row = start_line; row < end_line + 1; row++) {
      h_all_RowOffsets_[k] = SBM_permuted_csr->p[row] - SBM_permuted_csr->p[start_line];
      k++;
    }
  }

  // Copy this array to device
  int64_t *d_all_RowOffsets_ = NULL;
  CHECK_CUDA( cudaMalloc(&d_all_RowOffsets_, (2 * nc + n_colors) * sizeof(int64_t)) )
  CHECK_CUDA( cudaMemcpy(d_all_RowOffsets_, h_all_RowOffsets_, (2 * nc + n_colors) * sizeof(int64_t),
             cudaMemcpyHostToDevice) )

  int64_t* d_all_ColIndices_ = NULL;
  double* d_all_Values_ = NULL;
  size_t nnz = SBM_permuted_csr->p[2 * nc];
  printf("nnz = %d, nbblocks * 4 = %d\n", nnz, SBM_permuted->nbblocks * 4);

  CHECK_CUDA( cudaMalloc(&d_all_ColIndices_, nnz * sizeof(int64_t)) )
  CHECK_CUDA( cudaMalloc(&d_all_Values_, nnz * sizeof(double)) )

  CHECK_CUDA( cudaMemcpy(d_all_ColIndices_, SBM_permuted_csr->i, nnz * sizeof(int64_t),
             cudaMemcpyHostToDevice) )
  CHECK_CUDA( cudaMemcpy(d_all_Values_, SBM_permuted_csr->x, nnz * sizeof(double),
             cudaMemcpyHostToDevice) )

  cusparseConstSpMatDescr_t *all_cusparse_csrmat_ = (cusparseConstSpMatDescr_t *)malloc(n_colors * sizeof(cusparseConstSpMatDescr_t));
  
  size_t n_rows_;
  size_t current_nnz_;
  size_t cumul_nnz_ = 0;
  size_t cumul_n_rows_ = 0;
  for (unsigned int color = 0; color < n_colors; color++) {
    size_t start_line = 2 * sum_sizes[color];
    size_t end_line = 2 * sum_sizes[color + 1];
    n_rows_ = end_line - start_line;
    cumul_n_rows_ += n_rows_;
    current_nnz_ = SBM_permuted_csr->p[end_line] - SBM_permuted_csr->p[start_line];

    printf("%d\n", h_all_RowOffsets_[start_line + color]);

    CHECK_CUSPARSE( cusparseCreateConstCsr(&all_cusparse_csrmat_[color], n_rows_, 2 * nc, current_nnz_, 
                                           &d_all_RowOffsets_[start_line + color], &d_all_ColIndices_[cumul_nnz_], &d_all_Values_[cumul_nnz_], 
                                          CUSPARSE_INDEX_64I, CUSPARSE_INDEX_64I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) )

    cumul_nnz_ += current_nnz_;
  }

  printf("cumul_n_rows = %d, n_rows = %d\n", cumul_n_rows_, 2 * nc);
  printf("cumul_nnz = %d, nnz = %d\n", cumul_nnz_, SBM_permuted_csr->p[2 * sum_sizes[n_colors]]);

  printf("matrices created\n");
  fflush(stdout);

  // Create a res vector for each color
  double *h_all_res_ = (double *)calloc(2 * nc, sizeof(double));
  double *d_all_res_ = NULL;
  CHECK_CUDA( cudaMalloc(&d_all_res_, (2 * nc) * sizeof(double)) )
  CHECK_CUDA( cudaMemcpy(d_all_res_, h_all_res_, (2 * nc) * sizeof(double), cudaMemcpyHostToDevice) )

  cusparseDnVecDescr_t *all_vec_res_ = (cusparseDnVecDescr_t *)malloc(n_colors * sizeof(cusparseDnVecDescr_t));
  for (unsigned int color = 0; color < n_colors; color++) {
    CHECK_CUSPARSE( cusparseCreateDnVec(&all_vec_res_[color], 2 * (sum_sizes[color + 1] - sum_sizes[color]), &d_all_res_[2 * sum_sizes[color]], CUDA_R_64F) )
  }

  double* d_q = NULL;
  CHECK_CUDA( cudaMalloc(&d_q, (2 * nc) * sizeof(double)) )
  CHECK_CUDA( cudaMemcpy(d_q, problem->q, (2 * nc) * sizeof(double), cudaMemcpyHostToDevice) )
  cusparseDnVecDescr_t vecq;
  CHECK_CUSPARSE( cusparseCreateDnVec(&vecq, 2 * nc, d_q, CUDA_R_64F) )

  printf("vectors created\n");
  fflush(stdout);
  
  void **d_all_Buffer_ = (void **)malloc(n_colors * sizeof(void *));
  size_t *all_Buffer_Size_ = (size_t *)malloc(n_colors * sizeof(size_t));

  for (unsigned int color = 0; color < n_colors; color++) {
    CHECK_CUSPARSE( cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, all_cusparse_csrmat_[color],
                                          vecq, &beta, all_vec_res_[color], CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, &all_Buffer_Size_[color]) )
    CHECK_CUDA( cudaMalloc(&d_all_Buffer_[color], all_Buffer_Size_[color]) )

    CHECK_CUSPARSE( cusparseSpMV_preprocess(
                                handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha, all_cusparse_csrmat_[color], vecq, &beta, all_vec_res_[color], CUDA_R_64F,
                                CUSPARSE_SPMV_ALG_DEFAULT, d_all_Buffer_[color]) )

    CHECK_CUSPARSE( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, all_cusparse_csrmat_[color],
                                 vecq, &beta, all_vec_res_[color], CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, d_all_Buffer_[color]) )
  }
  

  CHECK_CUDA( cudaMemcpy(h_all_res_, d_all_res_, (2 * nc) * sizeof(double), cudaMemcpyDeviceToHost) )

  int correct = 1;
  for (int i = 0; i < 2 * nc; i++) {
      if (fabs(h_res_permuted[i] - h_all_res_[i]) / fabs(h_res_permuted[i]) > 1e-10) { // direct floating point comparison is not
        printf("%d %e %e %e\n", i, h_res_permuted[i], h_all_res_[i], fabs(h_res_permuted[i] - h_all_res_[i]) / fabs(h_res_permuted[i]));  
        correct = 0;             // reliable
          break;
      }
  }
  if (correct)
      printf("multicolor test PASSED\n");
  else
      printf("multicolor test FAILED\n"); */

  // Now initialize the right stuff

  // This only sets elements to 0, but I think they will still be present (as 0s) in the CSR conversion,
  // unless the conversion function checks for 0 elements
  remove_diagonal_blocks(problem); // set diagonal blocks of problem->M->matrix1 to 0

  CSparseMatrix *SBM_permuted_csr_nodiag = (CSparseMatrix *)malloc(sizeof(CSparseMatrix));
  size_t res = SBM_to_sparse_init_memory(problem->M->matrix1, SBM_permuted_csr_nodiag);
  res = SBM_to_sparse(problem->M->matrix1, SBM_permuted_csr_nodiag);
  if (res != 0) {
    printf("Error converting SBM to CSparse.\n");
    return;
  }
  else {
    printf("Conversion SBM->CSparse done\n");
  }

  // Create an array storing, for each color, the row offsets of the corresponding CSR submatrix 
  int64_t *h_all_RowOffsets = (int64_t *)malloc((2 * nc + n_colors) * sizeof(int64_t));
  size_t k = 0;
  for (unsigned int color = 0; color < n_colors; color++) {
    size_t start_line = 2 * sum_sizes[color];
    size_t end_line = 2 * sum_sizes[color + 1];

    for (unsigned int row = start_line; row < end_line + 1; row++) {
      h_all_RowOffsets[k] = SBM_permuted_csr_nodiag->p[row] - SBM_permuted_csr_nodiag->p[start_line];
      k++;
    }
  }

  // Copy this array to device
  int64_t *d_all_RowOffsets = NULL;
  CHECK_CUDA( cudaMalloc(&d_all_RowOffsets, (2 * nc + n_colors) * sizeof(int64_t)) )
  CHECK_CUDA( cudaMemcpy(d_all_RowOffsets, h_all_RowOffsets, (2 * nc + n_colors) * sizeof(int64_t), cudaMemcpyHostToDevice) )

  int64_t* d_all_ColIndices = NULL;
  double* d_all_Values = NULL;
  size_t nnz = SBM_permuted_csr_nodiag->p[2 * nc];

  CHECK_CUDA( cudaMalloc(&d_all_ColIndices, nnz * sizeof(int64_t)) )
  CHECK_CUDA( cudaMalloc(&d_all_Values, nnz * sizeof(double)) )

  CHECK_CUDA( cudaMemcpy(d_all_ColIndices, SBM_permuted_csr_nodiag->i, nnz * sizeof(int64_t),
             cudaMemcpyHostToDevice) )
  CHECK_CUDA( cudaMemcpy(d_all_Values, SBM_permuted_csr_nodiag->x, nnz * sizeof(double),
             cudaMemcpyHostToDevice) )

  cusparseConstSpMatDescr_t *all_cusparse_csrmat = (cusparseConstSpMatDescr_t *)malloc(n_colors * sizeof(cusparseConstSpMatDescr_t));
  
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
    current_nnz = SBM_permuted_csr_nodiag->p[end_line] - SBM_permuted_csr_nodiag->p[start_line];

    // printf("%d\n", h_all_RowOffsets[start_line + color]);

    CHECK_CUSPARSE( cusparseCreateConstCsr(&all_cusparse_csrmat[color], n_rows, 2 * nc, current_nnz, 
                                           &d_all_RowOffsets[start_line + color], &d_all_ColIndices[cumul_nnz], &d_all_Values[cumul_nnz], 
                                          CUSPARSE_INDEX_64I, CUSPARSE_INDEX_64I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) )

    cumul_nnz += current_nnz;
  }

  // Create vector for z, the reaction, which is the input of the product
  double *d_z = NULL;
  CHECK_CUDA( cudaMalloc(&d_z, (2 * nc) * sizeof(double)) )
  CHECK_CUDA( cudaMemcpy(d_z, z, (2 * nc) * sizeof(double), cudaMemcpyHostToDevice) )
  cusparseDnVecDescr_t vec_z;
  CHECK_CUSPARSE( cusparseCreateDnVec(&vec_z, 2 * nc, d_z, CUDA_R_64F) )
  // z vector created
  

  // Create the res vectors, which are the q's of the local problems
  double *h_q_new = (double *)malloc((2 * nc) * sizeof(double));
  double *d_q_new = NULL;
  CHECK_CUDA( cudaMalloc(&d_q_new, (2 * nc) * sizeof(double)) )
  CHECK_CUDA( cudaMemcpy(d_q_new, problem->q, (2 * nc) * sizeof(double), cudaMemcpyHostToDevice) )

  cusparseDnVecDescr_t *all_q_new = (cusparseDnVecDescr_t *)malloc(n_colors * sizeof(cusparseDnVecDescr_t));
  for (unsigned int color = 0; color < n_colors; color++) {
    CHECK_CUSPARSE( cusparseCreateDnVec(&all_q_new[color], 2 * (sum_sizes[color + 1] - sum_sizes[color]), &d_q_new[2 * sum_sizes[color]], CUDA_R_64F) )
  }
  
  void **d_all_Buffer = (void **)malloc(n_colors * sizeof(void *));
  size_t *all_Buffer_Size = (size_t *)malloc(n_colors * sizeof(size_t));

  alpha = 1.0;
  beta = 1.0;
  for (unsigned int color = 0; color < n_colors; color++) {
    CHECK_CUSPARSE( cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, all_cusparse_csrmat[color],
                                          vec_z, &beta, all_q_new[color], CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, &all_Buffer_Size[color]) )
    CHECK_CUDA( cudaMalloc(&d_all_Buffer[color], all_Buffer_Size[color]) )

    CHECK_CUSPARSE( cusparseSpMV_preprocess(
                                handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha, all_cusparse_csrmat[color], vec_z, &beta, all_q_new[color], CUDA_R_64F,
                                CUSPARSE_SPMV_ALG_DEFAULT, d_all_Buffer[color]) )
  }

  // Copy diagonal blocks to device
  double *d_diagonal_blocks = NULL;
  CHECK_CUDA( cudaMalloc(&d_diagonal_blocks, (4 * nc) * sizeof(double)) )
  CHECK_CUDA( cudaMemcpy(d_diagonal_blocks, diagonal_blocks, (4 * nc) * sizeof(double), cudaMemcpyHostToDevice) )

  // Copy determinants of diagonal blocks to device
  double *d_determinants = NULL;
  CHECK_CUDA( cudaMalloc(&d_determinants, nc * sizeof(double)) )
  CHECK_CUDA( cudaMemcpy(d_determinants, diagonal_block_determinant, nc * sizeof(double), cudaMemcpyHostToDevice) )

  // Copy mu to device
  double *d_mu = NULL;
  CHECK_CUDA( cudaMalloc(&d_mu, nc * sizeof(double)) )
  CHECK_CUDA( cudaMemcpy(d_mu, problem->mu, nc * sizeof(double), cudaMemcpyHostToDevice) )


  int threadsPerBlock = 256;
  int blocks = (nc + threadsPerBlock - 1) / threadsPerBlock;

  // REST 

  unsigned int nbblocks = problem->M->matrix1->nbblocks;
  double* blocks_contiguous = (double*)malloc(nbblocks * 4 * sizeof(double));
  double* current_block;
  for (unsigned int blockNum = 0; blockNum < nbblocks; blockNum++) {
    current_block = problem->M->matrix1->block[blockNum];
    for (unsigned int j = 0; j < 4; j++) {
      blocks_contiguous[blockNum * 4 + j] = current_block[j];
    }
  }

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

  size_t* index1_data;
  size_t* index2_data;

  unsigned int* freeze_contacts = NULL;
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
          for (unsigned int permuted_contact = sum_sizes[color];
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
    }  // end parallel region

    /***********************/
    /* NO FREEZING CONTACT */
    /***********************/
  } else {
    unsigned int contact;
    unsigned int pos;
    double light_error_sum = 0.;
    double localreaction[2];


    double* d_sumP2 = NULL;
    double* d_sumErr2 = NULL;

    CHECK_CUDA( cudaMalloc(&d_sumP2,   sizeof(double)) )
    CHECK_CUDA( cudaMalloc(&d_sumErr2, sizeof(double)) )

    CHECK_CUDA( cudaMemset(d_sumP2,   0, sizeof(double)) )
    CHECK_CUDA( cudaMemset(d_sumErr2, 0, sizeof(double)) )


    while ((iter < itermax) && has_not_converged) {

      // Set local problems q to q
      // maybe we can create a const q on device so that copy is faster?
      CHECK_CUDA( cudaMemcpy(d_q_new, problem->q, (2 * nc) * sizeof(double), cudaMemcpyHostToDevice) )
      CHECK_CUDA( cudaMemset(d_sumP2,   0, sizeof(double)) )
      CHECK_CUDA( cudaMemset(d_sumErr2, 0, sizeof(double)) )

      for (size_t color = 0; color < n_colors; color++) {
        int rangeSize = sum_sizes[color + 1] - sum_sizes[color];
        int threadsPerBlock = 256;
        int blocks = (rangeSize + threadsPerBlock - 1) / threadsPerBlock;
        size_t shmemSize = 2 * threadsPerBlock * sizeof(double);

        // Build local problems in parallel
        CHECK_CUSPARSE( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, all_cusparse_csrmat[color],
                                     vec_z, &beta, all_q_new[color], CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, d_all_Buffer[color]) )

        // Copy local problem q's back to host
        // Can I copy only the part that was modified? 
        // Like: CHECK_CUDA( cudaMemcpy(&h_q_new[...], &d_q_new[...], ... cudaMemcpyDeviceToHost) )
        fc2d_nsgs_local_solve_kernel_range_reduce<<<blocks, threadsPerBlock, shmemSize>>>(
          d_diagonal_blocks, d_determinants, d_q_new, d_mu, d_z,
          d_sumP2, d_sumErr2,
          sum_sizes[color], sum_sizes[color + 1]);
        // fc2d_nsgs_local_solve_kernel<<<blocks, threadsPerBlock>>>(d_diagonal_blocks, d_determinants, d_q_new, d_mu, d_z, sum_sizes[color], sum_sizes[color + 1]);


        /* CHECK_CUDA( cudaMemcpy(h_q_new, d_q_new, (2 * nc) * sizeof(double), cudaMemcpyDeviceToHost) )

        // Solve local problems (sequential for now)
        for (unsigned int permuted_contact = sum_sizes[color]; permuted_contact < sum_sizes[color + 1]; permuted_contact++) {
          pos = 2 * permuted_contact;
          localreaction[0] = z[pos];
          localreaction[1] = z[pos + 1];

          fc2d_nsgs_local_solve(&diagonal_blocks[4 * permuted_contact], diagonal_block_determinant[permuted_contact],
                                &h_q_new[pos], problem->mu[permuted_contact], localreaction);

          light_error_sum += light_error_squared(localreaction, &z[pos]);
          norm_r += squared_norm(localreaction);

          z[pos] = localreaction[0];
          z[pos + 1] = localreaction[1];
        } */

        CHECK_CUDA( cudaMemcpy(&light_error_sum, d_sumErr2, sizeof(double), cudaMemcpyDeviceToHost) )
        CHECK_CUDA( cudaMemcpy(&norm_r, d_sumP2, sizeof(double), cudaMemcpyDeviceToHost) )

        // Update z vector (same as above, can we only copy what was modified?)
        // CHECK_CUDA( cudaMemcpy(d_z, z, (2 * nc) * sizeof(double), cudaMemcpyHostToDevice) )

      }

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
  fc2d_free_diagonal_blocks(diagonal_blocks);
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

  free(SBM_col_permuted);
  free(SBM_permuted);
  free(q_permuted);
  free(mu_permuted);

  /* DO NOT FORGET TO FREE THE REST */
}
