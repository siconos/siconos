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

/*! \file SiconosAlgebraAddons.hpp
  \brief free functions for extra SiconosMatrix or Vectors operations
*/

#ifndef SICOALGEBRA_ADDONS_HPP
#define SICOALGEBRA_ADDONS_HPP

#include <concepts>
#include <optional>

#include "BlockVector.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace siconos::algebra {

// class BlockVector;

/** fill a CSparseMatrix (triplet) from non null values of a SiconosSparseMatrix
 *
 *  \param m input matrix
 *  \param csc the compressed column sparse matrix
 *  \param row_off
 *  \param col_off
 */
void fillTriplet(const SiconosSparseMatrix& m, CSparseMatrix* triplet, size_t row_off,
                 size_t col_off);

/** fill a CSparseMatrix (triplet) from non null values of a SiconosDenseMatrix
 *
 *  \param m input matrix
 *  \param csc the compressed column sparse matrix
 *  \param row_off
 *  \param col_off
 *  \param tol the tolerance to consider a number zero (not used if the matrix
 *  is sparse)
 */
void fillTriplet(SiconosDenseMatrix& m, CSparseMatrix* csc, size_t row_off, size_t col_off,
                 double tol = 1e-14);

/** \return a vector which contains the list of norm inf of each column of a matrix
    \param mat input matrix
 */
siconos::algebra::SiconosVector normInfByColumn(const SiconosDenseMatrix& m);

/** \return a vector which contains the list of norm inf of each column of a matrix
    \param mat input matrix
 */
siconos::algebra::SiconosVector normInfByColumn(const SiconosSparseMatrix& m);

/** \return true if the matrix is symmetric
    \param mat input matrix
    \param tol tolerance value
 */
bool isSymmetric(SiconosDenseMatrix& mat, double tol = 1e-12);
bool isSymmetric(const SiconosSparseMatrix& mat, double tol = 1e-12);

// Concept : any eigen matrix complient with middleCols * vec, which means for us dense, sparse
// and map(dense)
template <typename MatrixType>
concept MatrixBlockCompatible = requires(const MatrixType& A, Eigen::Index startCol,
                                         Eigen::Index blockSize, const SiconosVector& vec) {
  { A.middleCols(startCol, blockSize) * vec };  // -> std::same_as<SiconosVector>;
};

/**
 * Compute y += A.x with A a sparse or dense matrix and x a BlockVector
 *
 * @param[in] A input matrix (might be Dense, Map<Dense> or Sparse)
 * @param[in] x input vector (Block)
 * @param[in,out] y result
 * @param init true to reset y to zero before product, else accumulate in place
 */
template <MatrixBlockCompatible MatrixType>
void matrixBlockVector_prod(const MatrixType& A, const BlockVector& x,
                            Eigen::Ref<SiconosVector> y, bool init = true) {
  if (init) y.setZero();
  assert(y.size() == A.rows());
  assert(x.size() == A.cols());
  Eigen::Index startCol = 0;
  for (const auto& it : x) {
    Eigen::Index blockSize = it->size();
    auto subA = A.middleCols(startCol, blockSize);
    y.noalias() += subA * (*it);

    startCol += blockSize;
  }
}

// /** computes y = A*x or y += A*x if init = false
//   \param A a SiconosMatrix
//   \param x a SiconosVector
//   \param[in,out] y a SiconosVector
//   \param init a bool (default = true)
//   */
// void matrixVector_prod_toBlock(const SiconosMatrix& A, const SiconosVector& x, BlockVector&
// y,
//                                bool init = true);

/** computes y(block vector) = trans(A)*x (init = true) or y += trans(A)*x (init= false)
   works for dense and sparse
   \param x input vector
   \param A input matrix
   \param[in,out] y result
   \param init  false to accumulate result into y
  */
template <typename MatrixType>
void transposeMatrixVector_prod_toBlock(const SiconosVector& x, const MatrixType& A,
                                        BlockVector& y, bool init = true) {
  assert(A.rows() == x.size());
  assert(A.cols() == y.size());

  if (init) {
    for (auto& block : y) block->setZero();
  }

  Eigen::Index startCol = 0;
  for (auto& block : y) {
    Eigen::Index blockSize = block->size();
    auto subA = A.middleCols(startCol, blockSize);
    *block += subA.transpose() * x;
    startCol += blockSize;
  }
}

/** Generate a random sparse matrix. Useful for tests */
SiconosSparseMatrix generateRandomSparseMatrix(Eigen::Index rows, Eigen::Index cols, int nnz,
                                               std::optional<double> density = std::nullopt);
}  // namespace siconos::algebra

#endif
