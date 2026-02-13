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
#include <utility>

siconos::algebra::SiconosVector siconos::algebra::normInfByColumn(
    const siconos::algebra::SiconosDenseMatrix& m) {
  siconos::algebra::SiconosVector v(m.cols());
  for (siconos::algebra::Index j = 0; j < m.cols(); ++j)
    v(j) = m.col(j).lpNorm<Eigen::Infinity>();
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

bool siconos::algebra::isSymmetric(const SiconosDenseMatrix& mat, double tol) {
  // works with sparse and dense
  // Warning: might be quite expensive
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

  const SparseIndex* outer = m.outerIndexPtr();  // column pointers
  const SparseIndex* inner = m.innerIndexPtr();  // row indices
  const double* values = m.valuePtr();           // non-zero values

  Index cols = m.cols();

  for (Index j = 0; j < cols; ++j) {
    for (auto idx = outer[j]; idx < outer[j + 1]; ++idx) {
      // Insert triplet directly, no threshold
      size_t idx_s = static_cast<size_t>(idx);
      cs_entry(triplet, inner[idx_s] + row_off, j + col_off, values[idx_s]);
    }
  }
}

void siconos::algebra::fillTriplet(SiconosDenseMatrix& m, CSparseMatrix* triplet,
                                   size_t row_off, size_t col_off, double tol) {
  assert(triplet);
  size_t nrow = static_cast<size_t>(m.rows());
  size_t ncol = static_cast<size_t>(m.cols());

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
    siconos::algebra::Index rows, siconos::algebra::Index cols, SparseIndex nnz,
    std::optional<double> density) {
  std::mt19937 rng(42);
  std::uniform_int_distribution<siconos::algebra::Index> rowDist(0, rows - 1);
  std::uniform_int_distribution<siconos::algebra::Index> colDist(0, cols - 1);
  std::uniform_real_distribution<double> valDist(1e-6, 1.0);  // never equal to zero
  std::bernoulli_distribution signDist(0.5);

  if (density) {
    assert(*density >= 0.0 && *density <= 1.0);
    nnz = static_cast<SparseIndex>(*density * rows) * cols;
  }

  // --- Generate keys (i,j) encoded in an uint64_t
  std::vector<uint64_t> keys;
  assert(nnz >= 0);
  keys.reserve(static_cast<size_t>(nnz * 1.1));
  while (std::ssize(keys) < nnz) {
    uint64_t key =
        (static_cast<uint64_t>(rowDist(rng)) << 32) | static_cast<uint32_t>(colDist(rng));
    keys.push_back(key);
  }

  // --- Sort and eliminate duplicate val
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());

  // --- Add values to match nnz
  while (std::ssize(keys) < nnz) {
    uint64_t key =
        (static_cast<uint64_t>(rowDist(rng)) << 32) | static_cast<uint32_t>(colDist(rng));
    if (!std::binary_search(keys.begin(), keys.end(), key)) {
      keys.push_back(key);
      std::inplace_merge(keys.begin(), keys.end() - 1, keys.end());  // to keep the same order
    }
  }

  // --- Generate triplets
  std::vector<Eigen::Triplet<double>> triplets;
  assert(nnz >= 0);
  triplets.reserve(static_cast<size_t>(nnz));

  for (uint64_t key : keys) {
    Index i = static_cast<Index>(key >> 32);
    Index j = static_cast<Index>(key & 0xffffffffu);
    double val = (signDist(rng) ? 1.0 : -1.0) * valDist(rng);
    triplets.emplace_back(i, j, val);
  }
  siconos::algebra::SiconosSparseMatrix mat(rows, cols);
  mat.setFromTriplets(triplets.begin(), triplets.end());
  return mat;
}