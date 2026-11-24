/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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

#include "plasticity_2d_stress_coordinates.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "NumericsMatrix.h"
#include "SparseBlockMatrix.h"

#define H_ M_SQRT1_2

/* The four 3x3 maps (row-major), for Mandel and Voigt notations:
 *   P (s -> sigma_c), P^{-1} (sigma_c -> s), P^{-T} (w_s -> w_c), P^T (w_c -> w_s).
 * Mandel: P orthogonal, P^{-1} = P^T, P^{-T} = P.  Voigt: P = P_M D, D = diag(1, 1, sqrt2). */
enum { MAP_P = 0, MAP_P_INV = 1, MAP_P_INV_T = 2, MAP_P_T = 3 };

static const double MANDEL[4][9] = {
    {-H_, -H_, 0., H_, -H_, 0., 0., 0., 1.}, /* P */
    {-H_, H_, 0., -H_, -H_, 0., 0., 0., 1.}, /* P^{-1} = P^T */
    {-H_, -H_, 0., H_, -H_, 0., 0., 0., 1.}, /* P^{-T} = P */
    {-H_, H_, 0., -H_, -H_, 0., 0., 0., 1.}  /* P^T */
};
static const double VOIGT[4][9] = {
    {-H_, -H_, 0., H_, -H_, 0., 0., 0., M_SQRT2}, /* P = P_M D */
    {-H_, H_, 0., -H_, -H_, 0., 0., 0., H_},      /* P^{-1} = D^{-1} P_M^T */
    {-H_, -H_, 0., H_, -H_, 0., 0., 0., H_},      /* P^{-T} = P_M D^{-1} */
    {-H_, H_, 0., -H_, -H_, 0., 0., 0., M_SQRT2}  /* P^T = D P_M^T */
};

static const double* get_map(PlasticityStressCoordinates coordinates, int which) {
  switch (coordinates) {
    case PLASTICITY_STRESS_COORDINATES_CARTESIAN_MANDEL:
      return MANDEL[which];
    case PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT:
      return VOIGT[which];
    default:
      return NULL; /* cone coordinates: identity */
  }
}

static void apply_map(const double* A, int nc, const double* x, double* y) {
  if (!A) {
    if (x != y) memcpy(y, x, (size_t)(3 * nc) * sizeof(double));
    return;
  }
  for (int e = 0; e < nc; e++) {
    const double* xe = x + 3 * e;
    double t0 = A[0] * xe[0] + A[1] * xe[1] + A[2] * xe[2];
    double t1 = A[3] * xe[0] + A[4] * xe[1] + A[5] * xe[2];
    double t2 = A[6] * xe[0] + A[7] * xe[1] + A[8] * xe[2];
    double* ye = y + 3 * e;
    ye[0] = t0;
    ye[1] = t1;
    ye[2] = t2;
  }
}

void plasticity_2d_stress_to_cone(PlasticityStressCoordinates coordinates, int nc,
                                  const double* s, double* sigma_c) {
  apply_map(get_map(coordinates, MAP_P), nc, s, sigma_c);
}

void plasticity_2d_stress_from_cone(PlasticityStressCoordinates coordinates, int nc,
                                    const double* sigma_c, double* s) {
  apply_map(get_map(coordinates, MAP_P_INV), nc, sigma_c, s);
}

void plasticity_2d_strainrate_to_cone(PlasticityStressCoordinates coordinates, int nc,
                                      const double* w_s, double* w_c) {
  apply_map(get_map(coordinates, MAP_P_INV_T), nc, w_s, w_c);
}

void plasticity_2d_strainrate_from_cone(PlasticityStressCoordinates coordinates, int nc,
                                        const double* w_c, double* w_s) {
  apply_map(get_map(coordinates, MAP_P_T), nc, w_c, w_s);
}

/* X <- A X B for the 3x3 column-major sub-matrix X (leading dimension ld), A and B
 * row-major. */
static void transform_block3(const double* A, const double* B, double* X, size_t ld) {
  double AX[9]; /* row-major */
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      AX[3 * i + j] = A[3 * i] * X[0 + j * ld] + A[3 * i + 1] * X[1 + j * ld] +
                      A[3 * i + 2] * X[2 + j * ld];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      X[i + j * ld] = AX[3 * i] * B[j] + AX[3 * i + 1] * B[3 + j] + AX[3 * i + 2] * B[6 + j];
}

/* X <- A X B on all the 3x3 sub-blocks of a column-major (rows x cols) matrix. */
static void transform_dense(const double* A, const double* B, double* X, size_t rows,
                            size_t cols) {
  for (size_t j = 0; j < cols; j += 3)
    for (size_t i = 0; i < rows; i += 3) transform_block3(A, B, X + i + j * rows, rows);
}

/* block diagonal matrix with the 3x3 block A (row-major), NM_SPARSE */
static NumericsMatrix* block_diagonal(const double* A, int n) {
  NumericsMatrix* T = NM_create(NM_SPARSE, n, n);
  NM_triplet_alloc(T, 9 * (n / 3));
  for (int e = 0; e < n; e += 3)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        if (A[3 * i + j] != 0.0) NM_entry(T, (size_t)(e + i), (size_t)(e + j), A[3 * i + j]);
  return T;
}

static int sbm_blocks_multiple_of_3(const SparseBlockStructuredMatrix* B) {
  for (size_t i = 0; i < B->blocknumber0; i++)
    if (B->blocksize0[i] % 3) return 0;
  for (size_t j = 0; j < B->blocknumber1; j++)
    if (B->blocksize1[j] % 3) return 0;
  return 1;
}

NumericsMatrix* plasticity_2d_matrix_to_cone(PlasticityStressCoordinates coordinates,
                                             NumericsMatrix* M) {
  const double* A = get_map(coordinates, MAP_P_INV_T);
  const double* B = get_map(coordinates, MAP_P_INV);
  if (!A) return NM_copy(M);

  if (M->storageType == NM_DENSE) {
    NumericsMatrix* Mc = NM_copy(M);
    transform_dense(A, B, Mc->dense, (size_t)Mc->size0, (size_t)Mc->size1);
    return Mc;
  }
  if (M->storageType == NM_SPARSE_BLOCK && sbm_blocks_multiple_of_3(M->sparse_block)) {
    NumericsMatrix* Mc = NM_copy(M);
    SparseBlockStructuredMatrix* S = Mc->sparse_block;
    /* blocksize0/1 are cumulated sizes */
    for (size_t bi = 0; bi < S->filled1 - 1; bi++) {
      size_t rows = S->blocksize0[bi] - (bi ? S->blocksize0[bi - 1] : 0);
      for (size_t k = S->index1_data[bi]; k < S->index1_data[bi + 1]; k++) {
        size_t bj = S->index2_data[k];
        size_t cols = S->blocksize1[bj] - (bj ? S->blocksize1[bj - 1] : 0);
        transform_dense(A, B, S->block[k], rows, cols);
      }
    }
    return Mc;
  }
  /* any other storage: products with sparse block-diagonal matrices */
  int n = M->size0;
  NumericsMatrix* TA = block_diagonal(A, n);
  NumericsMatrix* TB = block_diagonal(B, n);
  NumericsMatrix* MB = NM_multiply(M, TB);
  NumericsMatrix* Mc = NM_multiply(TA, MB);
  NM_free(MB);
  NM_free(TA);
  NM_free(TB);
  return Mc;
}

PlasticityProblem* plasticity_2d_problem_to_cone(PlasticityProblem* problem) {
  PlasticityProblem* cone = plasticity_copy(problem);
  if (problem->stress_coordinates == PLASTICITY_STRESS_COORDINATES_CONE) return cone;
  NM_free(cone->M);
  cone->M = plasticity_2d_matrix_to_cone(problem->stress_coordinates, problem->M);
  plasticity_2d_strainrate_to_cone(problem->stress_coordinates, problem->numberOfCones,
                                   problem->q, cone->q);
  cone->stress_coordinates = PLASTICITY_STRESS_COORDINATES_CONE;
  return cone;
}
