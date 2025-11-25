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
#include "fc2d_cuda_kernels.h"
#include <math.h>

/**
 * Device function: Solve 2x2 friction contact local problem.
 * Based on fc2d_nsgs_local_solve from fc2d_nsgs_graph_permut.c (lines 241-285).
 */
__device__ void fc2d_local_solve_device(const double *W, double D, const double *q, double mu,
                                        double *P) {
  /* Unpack W matrix elements (column-major storage) */
  double Wnn = W[0];
  double Wtn = W[1];
  double Wnt = W[2];
  double Wtt = W[3];

  if(q[0] > 0.0) {
    /* No contact */
    P[0] = 0.0;
    P[1] = 0.0;
  } else {
    /* Solve WP + q = 0 */
    P[0] = -(Wtt * q[0] - Wnt * q[1]) / D;
    P[1] = -(-Wtn * q[0] + Wnn * q[1]) / D;

    double muPn = mu * P[0];

    if(fabs(P[1]) > muPn) {
      /* Outside friction cone */
      if(P[1] + muPn < 0.0) {
        P[0] = -q[0] / (Wnn - mu * Wnt);
        P[1] = -mu * P[0];
      } else {
        P[0] = -q[0] / (Wnn + mu * Wnt);
        P[1] = mu * P[0];
      }
    }
  }
}

/**
 * Device function: Build local problem for contact i.
 * Computes q_local = q + M_{i,:} * z (excluding diagonal block).
 * Based on fc2d_nsgs_buildLocalProblem_parallel from fc2d_nsgs_graph_permut.c (lines 73-101).
 */
__device__ void fc2d_build_local_problem_device(unsigned int contact, const int *bsrRowPtr,
                                                const int *bsrColInd, const double *bsrVal,
                                                const double *z, const double *q_base,
                                                double *q_local) {
  /* Initialize with base q */
  q_local[0] = q_base[2 * contact];
  q_local[1] = q_base[2 * contact + 1];

  /* Add M_{i,:} * z (excluding diagonal) */
  int row_start = bsrRowPtr[contact];
  int row_end = bsrRowPtr[contact + 1];

  for(int blockNum = row_start; blockNum < row_end; blockNum++) {
    int colNum = bsrColInd[blockNum];

    if(colNum != (int)contact) {
      /* 2x2 block matrix-vector product (block stored column-major) */
      const double *block = &bsrVal[blockNum * 4];
      const double *z_col = &z[2 * colNum];

      /* q_local += block * z_col */
      q_local[0] += block[0] * z_col[0] + block[2] * z_col[1];
      q_local[1] += block[1] * z_col[0] + block[3] * z_col[1];
    }
  }
}

/**
 * Main GPU kernel: Process one color of contacts in Gauss-Seidel iteration.
 * Each thread handles one contact.
 */
__global__ void fc2d_nsgs_color_iteration_kernel(unsigned int contact_start,
                                                  unsigned int contact_end, const int *bsrRowPtr,
                                                  const int *bsrColInd, const double *bsrVal,
                                                  const double *diagonal_blocks,
                                                  const double *diagonal_determinants,
                                                  const double *q_base, const double *mu,
                                                  double *z, double *error_local,
                                                  double *norm_local) {
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int contact = contact_start + idx;

  if(contact >= contact_end)
    return;

  /* Local variables */
  double q_local[2];
  double z_old[2];
  double z_new[2];

  /* Save old z */
  unsigned int pos = 2 * contact;
  z_old[0] = z[pos];
  z_old[1] = z[pos + 1];

  /* Build local problem */
  fc2d_build_local_problem_device(contact, bsrRowPtr, bsrColInd, bsrVal, z, q_base, q_local);

  /* Solve local problem */
  const double *W = &diagonal_blocks[4 * contact];
  double D = diagonal_determinants[contact];
  fc2d_local_solve_device(W, D, q_local, mu[contact], z_new);

  /* Update z */
  z[pos] = z_new[0];
  z[pos + 1] = z_new[1];

  /* Compute error contribution */
  double dx0 = z_new[0] - z_old[0];
  double dx1 = z_new[1] - z_old[1];
  error_local[idx] = dx0 * dx0 + dx1 * dx1;
  norm_local[idx] = z_new[0] * z_new[0] + z_new[1] * z_new[1];
}

/**
 * Host function: Launch kernel for one color iteration.
 */
void fc2d_nsgs_color_iteration(unsigned int contact_start, unsigned int contact_end,
                                const int *d_bsrRowPtr, const int *d_bsrColInd,
                                const double *d_bsrVal, const double *d_diagonal_blocks,
                                const double *d_diagonal_determinants, const double *d_q,
                                const double *d_mu, double *d_z, double *d_error_local,
                                double *d_norm_local, cudaStream_t stream) {
  unsigned int num_contacts = contact_end - contact_start;

  /* Launch configuration */
  const int threads_per_block = 256;
  int num_blocks = (num_contacts + threads_per_block - 1) / threads_per_block;

  /* Launch kernel */
  fc2d_nsgs_color_iteration_kernel<<<num_blocks, threads_per_block, 0, stream>>>(
      contact_start, contact_end, d_bsrRowPtr, d_bsrColInd, d_bsrVal, d_diagonal_blocks,
      d_diagonal_determinants, d_q, d_mu, d_z, d_error_local, d_norm_local);

  CUDA_CHECK(cudaGetLastError());
}

/**
 * GPU kernel: Reduction for sum computation using shared memory.
 */
__global__ void reduce_sum_kernel(const double *input, double *output, unsigned int n) {
  extern __shared__ double sdata[];

  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  /* Load data into shared memory */
  sdata[tid] = (i < n) ? input[i] : 0.0;
  __syncthreads();

  /* Reduction in shared memory */
  for(unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if(tid < s) {
      sdata[tid] += sdata[tid + s];
    }
    __syncthreads();
  }

  /* Write result for this block using atomic add */
  if(tid == 0) {
    atomicAdd(output, sdata[0]);
  }
}

/**
 * Host function: Reduce error and norm arrays and copy to host.
 */
void reduce_error_norm(const double *d_error_local, const double *d_norm_local, unsigned int n,
                       double *h_error, double *h_norm) {
  /* Allocate device output */
  double *d_error_sum, *d_norm_sum;
  CUDA_CHECK(cudaMalloc(&d_error_sum, sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_norm_sum, sizeof(double)));
  CUDA_CHECK(cudaMemset(d_error_sum, 0, sizeof(double)));
  CUDA_CHECK(cudaMemset(d_norm_sum, 0, sizeof(double)));

  /* Launch reduction */
  const int threads = 256;
  int blocks = (n + threads - 1) / threads;
  size_t smem_size = threads * sizeof(double);

  reduce_sum_kernel<<<blocks, threads, smem_size>>>(d_error_local, d_error_sum, n);
  reduce_sum_kernel<<<blocks, threads, smem_size>>>(d_norm_local, d_norm_sum, n);

  /* Copy results to host */
  CUDA_CHECK(cudaMemcpy(h_error, d_error_sum, sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(h_norm, d_norm_sum, sizeof(double), cudaMemcpyDeviceToHost));

  /* Cleanup */
  CUDA_CHECK(cudaFree(d_error_sum));
  CUDA_CHECK(cudaFree(d_norm_sum));
}
