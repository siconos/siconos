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

// void siconos::algebra::matrixVector_prod_toBlock(const SiconosMatrix& A,
//                                                  const SiconosVector& x, BlockVector& y,
//                                                  bool init) {
//   unsigned int startRow = 0;
//   // For Each subvector of y, y[i], private_prod computes y[i] = subA x, subA
//   // being a submatrix of A corresponding to y[i] position.
//   //       // private_prod takes into account the fact that x and y[i] may be
//   if (init) {
//     for (auto& it : y) {
//       it->setZero();
//     }
//   }
//   //       block vectors.
//   for (auto& it : y) {
//     it->noalias() += A.block(startRow, 0, it->size(), x.size()) * x;
//     startRow += it->size();
//   }
// }

void siconos::algebra::matrixBlockVector_prod(const SiconosMatrix& A, const BlockVector& x,
                                              SiconosVector& y, bool init) {
  if (init) y.setZero();
  unsigned int startRow = 0;
  unsigned int startCol = 0;
  // In private_addprod, the sum of all blocks of x, x[i], is computed: y =
  // Sum_i (subA x[i]), with subA a submatrix of A, starting from position
  // startRow in rows and startCol in columns. private_prod takes also into
  // account the fact that each block of x can also be a block.
  for (auto& it : x) {
    assert(&y != &*it);
    y += A.block(startRow, startCol, y.size(), it->size()) * *it;
    startCol += it->size();
  }
}

void siconos::algebra::matrixBlockVector_prod(const SiconosMatrix& A, const BlockVector& x,
                                              Eigen::Ref<SiconosVector> y, bool init) {
  if (init) y.setZero();
  unsigned int startRow = 0;
  unsigned int startCol = 0;
  // In private_addprod, the sum of all blocks of x, x[i], is computed: y =
  // Sum_i (subA x[i]), with subA a submatrix of A, starting from position
  // startRow in rows and startCol in columns. private_prod takes also into
  // account the fact that each block of x can also be a block.
  for (auto& it : x) {
    y += A.block(startRow, startCol, y.size(), it->size()) * *it;
    startCol += it->size();
  }
}

void siconos::algebra::transposeMatrixVector_prod_toBlock(const SiconosVector& x,
                                                          const SiconosMatrix& A,
                                                          BlockVector& y, bool init) {
  if (A.rows() != x.size()) THROW_EXCEPTION("inconsistent sizes between A and x.");

  if (A.cols() != y.size()) THROW_EXCEPTION("inconsistent sizes between A and y.");
  if (init) {  // y = subA * x , else y += subA * x
    for (auto& it : y) it->setZero();
  }

  unsigned int pos = 0;
  // For Each subvector of y, y[i], computes y[i] = transpose(subA) x, subA
  // being a submatrix of A corresponding to y[i] position. private_prod takes
  // into account the fact that x and y[i] may be block vectors.
  auto sizeX = x.size();
  for (auto& it : y) {
    // we take a submatrix subA of A, starting from row startRow to row
    // (startRow+sizeY) and between columns startCol and (startCol+sizeX). Then
    // computation of y = subA*x + y.

    auto sizeY = it->size();

    assert(&*it != &x);
    it->noalias() += A.transpose().block(pos, 0, sizeY, sizeX) * x;
    pos += it->size();
  }
}