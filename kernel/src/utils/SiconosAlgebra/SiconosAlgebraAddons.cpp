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

#include "SiconosAlgebraAddons.hpp"

#include <random>
#include <set>

#include "BlockVector.hpp"
#include "SiconosException.hpp"
#include "Tools.hpp"  // toString

siconos::algebra::SiconosVector siconos::algebra::normInfByColumn(
    const siconos::algebra::SiconosDenseMatrix& m) {
  siconos::algebra::SiconosVector v(m.cols());
  for (Eigen::Index j = 0; j < m.cols(); ++j) v(j) = m.col(j).lpNorm<Eigen::Infinity>();
  return v;
}

siconos::algebra::SiconosVector siconos::algebra::normInfByColumn(
    const siconos::algebra::SiconosSparseMatrix& m) {
  SiconosVector v = SiconosVector::Zero(m.cols());
  for (int k = 0; k < m.outerSize(); ++k) {
    for (SiconosSparseMatrix::InnerIterator it(m, k); it; ++it) {
      double absVal = std::abs(it.value());
      if (absVal > v(it.col())) v(it.col()) = absVal;
    }
  }
  return v;  // RVO
}

bool siconos::algebra::isSymmetric(SiconosDenseMatrix& mat, double tol) {
  // works with sparse and dense
  // Warning: might be quite expensive, especially for sparse matrices
  if (mat.rows() != mat.cols()) return false;
  return mat.isApprox(mat.transpose(), tol);
}

bool siconos::algebra::isSymmetric(const SiconosSparseMatrix& mat, double tol) {
  if (mat.rows() != mat.cols()) return false;
  for (int col = 0; col < mat.outerSize(); ++col) {
    for (SiconosSparseMatrix::InnerIterator it(mat, col); it; ++it) {
      int row = it.row();
      double a_ij = it.value();
      double a_ji = mat.coeff(col, row);  // A(j,i)
      if (std::abs(a_ij - a_ji) > tol) return false;
    }
  }
  return true;
}

void siconos::algebra::fillTriplet(const SiconosSparseMatrix& m, CSparseMatrix* triplet,
                                   size_t row_off, size_t col_off) {
  assert(triplet);

  const int* outer = m.outerIndexPtr();  // column pointers
  const int* inner = m.innerIndexPtr();  // row indices
  const double* values = m.valuePtr();   // non-zero values

  int cols = m.cols();

  for (int j = 0; j < cols; ++j) {
    for (int idx = outer[j]; idx < outer[j + 1]; ++idx) {
      // Insert triplet directly, no threshold check
      cs_entry(triplet, inner[idx] + row_off, j + col_off, values[idx]);
    }
  }
}

void siconos::algebra::fillTriplet(SiconosDenseMatrix& m, CSparseMatrix* triplet,
                                   size_t row_off, size_t col_off, double tol) {
  assert(triplet);
  size_t nrow = m.rows();
  size_t ncol = m.cols();

  double* arr = m.data();
  for (size_t j = 0; j < ncol; ++j) {
    for (size_t i = 0; i < nrow; ++i) {
      // col-major

      CSparseMatrix_zentry(triplet, i + row_off, j + col_off, arr[i + j * nrow],
                           std::numeric_limits<double>::epsilon());
    }
  }
}

siconos::algebra::SiconosSparseMatrix siconos::algebra::generateRandomSparseMatrix(
    Eigen::Index rows, Eigen::Index cols, int nnz, std::optional<double> density) {
  siconos::algebra::SiconosSparseMatrix mat(rows, cols);
  std::vector<Eigen::Triplet<double>> triplets;

  std::mt19937 rng(42);
  std::uniform_int_distribution<int> rowDist(0, rows - 1);
  std::uniform_int_distribution<int> colDist(0, cols - 1);
  std::uniform_real_distribution<double> valDist(1e-6, 1.0);  // never equal to zero
  std::bernoulli_distribution signDist(0.5);

  if (density) {
    assert(*density <= 1.);
    assert(*density >= 0.);
    nnz = static_cast<int>(rows * cols * *density);
  }

  triplets.reserve(nnz);
  std::set<std::pair<int, int>> used_positions;

  // Step 1: Generate unique positions
  while ((int)used_positions.size() < nnz) {
    int i = rowDist(rng);
    int j = colDist(rng);
    used_positions.emplace(i, j);
  }

  // Step 2: Generate triplets with non-zero values
  for (const auto& [i, j] : used_positions) {
    double val = (signDist(rng) ? 1.0 : -1.0) * valDist(rng);
    triplets.emplace_back(i, j, val);
  }

  mat.setFromTriplets(triplets.begin(), triplets.end());
  return mat;
}