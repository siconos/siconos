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

#ifndef FC2D_CUDA_KERNELS_H
#define FC2D_CUDA_KERNELS_H

#include <cuda_runtime.h>

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/**
 * Launch kernel to process one color of contacts in the Gauss-Seidel iteration.
 *
 * \param contact_start First contact index in this color
 * \param contact_end Last contact index (exclusive) in this color
 * \param d_bsrRowPtr Device pointer to BSR row pointers
 * \param d_bsrColInd Device pointer to BSR column indices
 * \param d_bsrVal Device pointer to BSR block values
 * \param d_diagonal_blocks Device pointer to diagonal blocks (2x2)
 * \param d_diagonal_determinants Device pointer to precomputed determinants
 * \param d_q Device pointer to right-hand side vector
 * \param d_mu Device pointer to friction coefficients
 * \param d_z Device pointer to reaction vector (input/output)
 * \param d_error_local Device pointer to per-thread error contributions (output)
 * \param d_norm_local Device pointer to per-thread norm contributions (output)
 * \param stream CUDA stream for asynchronous execution
 */
void fc2d_nsgs_color_iteration(unsigned int contact_start, unsigned int contact_end,
                                const int *d_bsrRowPtr, const int *d_bsrColInd,
                                const double *d_bsrVal, const double *d_diagonal_blocks,
                                const double *d_diagonal_determinants, const double *d_q,
                                const double *d_mu, double *d_z, double *d_error_local,
                                double *d_norm_local, cudaStream_t stream);

/**
 * Reduce error and norm arrays on GPU and copy results to host.
 *
 * \param d_error_local Device pointer to per-thread error contributions
 * \param d_norm_local Device pointer to per-thread norm contributions
 * \param n Number of elements to reduce
 * \param h_error Host pointer to store reduced error sum
 * \param h_norm Host pointer to store reduced norm sum
 */
void reduce_error_norm(const double *d_error_local, const double *d_norm_local, unsigned int n,
                       double *h_error, double *h_norm);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif /* FC2D_CUDA_KERNELS_H */
