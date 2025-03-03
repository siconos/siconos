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

/*! \file SiconosMatrix.hpp
  Interface for matrices handling.
*/

#ifndef SICOMAT
#define SICOMAT

#include "EigenInclude.hpp"  // Must be included before Eigen/Core
//
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/SparseCore>  // For Eigen Sparse matrices

#include "CSparseMatrix.h"  // For CSparseMatrix
#include "SiconosVector.hpp"

struct NumericsMatrix;

namespace siconos::algebra {

/**
   Abstract class to provide interface for matrices handling

   Matrices can be either block or Simple.
   See Derived classes for details.

   In Siconos, a "matrix" can be either a SiconosMatrix or a BlockMatrix, ie a
   container of several pointers to SiconosMatrix

   You can find an overview on how to build and use vectors and matrices in
   siconos users guide
   .

*/
using SiconosDenseMatrix =
    Eigen::Matrix<double_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using SiconosMatrix33 = Eigen::Matrix<double_t, 3, 3, Eigen::ColMajor>;
/** Sparse matrix storage */
using SiconosSparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
using Triplet = Eigen::Triplet<double>;  // Used to fill sparse matrices

// Select between dense and sparse storage
// SICONOS_SPARSE is set by cmake.
// Use -DSICONOS_SPARSE=1 to activate sparse storage
#ifdef SICONOS_SPARSE
using SiconosMatrix = SiconosSparseMatrix;
#else
using SiconosMatrix = SiconosDenseMatrix;
#endif
// Map types

using MapVectorType = Eigen::Map<SiconosVector>;
using MapVector3Type = Eigen::Map<SiconosVector3>;
using ConstMapVectorType = Eigen::Map<const SiconosVector>;
using ConstMapVector3Type = Eigen::Map<const SiconosVector3>;
using MapType = Eigen::Map<SiconosMatrix>;
using ConstMapType = Eigen::Map<const SiconosMatrix>;

using SiconosDiagonalMatrix = Eigen::DiagonalMatrix<double_t, Eigen::Dynamic>;

// Tmp version. To replace BlockVector ?
using EigenBlock = std::vector<Eigen::Ref<SiconosVector>>;

// type used as return value for vectors and matrix dimensions.
using SiconosSize_t = Eigen::Index;

enum class StorageType { dense, sparse };

/** return the number of non-zero in the matrix
 *
 *  \param csc the compressed column sparse matrix
 *  \param row_off
 *  \param col_off
 *  \param tol the tolerance to consider a number zero (not used if the matrix
 * is sparse) \return the number of non-zeros
 */
bool fillTriplet(SiconosMatrix &m, CSparseMatrix *csc, size_t row_off, size_t col_off,
                 double tol = 1e-14);
void normInfByColumn(const SiconosMatrix &m, SiconosVector &v);

bool checkSymmetry(SiconosMatrix &m, double tol);

siconos::algebra::SiconosMatrix readMatrixFromFile(const std::string &filename,
                                                   bool ascii = true);

/** General dense matrix print function
 *  \param mat the matrix to be displayed
 */
template <typename Derived>
void print(const Eigen::EigenBase<Derived> &mat) {
  std::cout << std::scientific << std::setprecision(6) << mat.derived() << std::endl;
}

// Specialization for sparse matrices
inline void print(const SiconosSparseMatrix &mat) {
  std::cout << "Sparse Matrix (" << mat.rows() << "x" << mat.cols()
            << "), non-zeros: " << mat.nonZeros() << "\n";

  for (int k = 0; k < mat.outerSize(); ++k) {
    for (SiconosSparseMatrix::InnerIterator it(mat, k); it; ++it) {
      std::cout << "(" << it.row() << ", " << it.col() << ") -> " << std::scientific
                << std::setprecision(6) << it.value() << "\n";
    }
  }
}

// Specialization for column vectors (prints as row)
inline void print(const SiconosVector &vec) {
  std::cout << std::scientific << std::setprecision(6) << vec.transpose() << "\n";
}

/** \return inf. norm of the matrix
 *  \param mat the input matrix
 */
template <typename Derived>
auto normInf(const Eigen::MatrixBase<Derived> &mat) {
  return mat.cwiseAbs().rowwise().sum().maxCoeff();
}

// Specialization for sparse matrices (uses iterators)
inline double normInf(const SiconosSparseMatrix &mat) {
  double maxNorm = 0.0;
  // Warn: col-major
  std::vector<double> rowSums(mat.rows(), 0.0);

  for (int k = 0; k < mat.outerSize(); ++k) {
    // Sum over each row
    for (SiconosSparseMatrix::InnerIterator it(mat, k); it; ++it) {
      rowSums[it.row()] += std::abs(it.value());
    }
  }
  for (double rowSum : rowSums) {
    maxNorm = std::max(maxNorm, rowSum);
  }
  return maxNorm;
}
template <typename Derived>
void normInf(const Eigen::EigenBase<Derived> &mat) {
  return mat->cwiseAbs().rowwise().sum().maxCoeff();
}

/** Set a sub-block of a matrix
 *  \param row_min starting row index of the written block
 *  \param col_min starting col index of the written block
 *  \param[in] block_input input matrix to be written in the block
 *  \param[in,out] m_out matrix to be updated
 */
template <typename Derived>
inline void setBlock(typename Derived::Index row_min, typename Derived::Index col_min,
                     const Eigen::MatrixBase<Derived> &block_input,
                     Eigen::Ref<SiconosMatrix> m_out) {
  // assert(m != *this);
  assert(row_min < m_out.rows() && "row is out of range");
  assert(col_min < m_out.cols() && "column is out of range");

  // Get the number of rows and columns in the input block
  auto block_rows = block_input.rows();
  auto block_cols = block_input.cols();

  // Calculate the max row and column indices
  assert(row_min + block_rows <= m_out.rows() && "block goes out of range in rows");
  assert(col_min + block_cols <= m_out.cols() && "block goes out of range in columns");

  m_out.block(row_min, col_min, block_rows, block_cols) = block_input;
}

/** Set a sub-block of a Vector
 *  \param row_min starting row index of the written block
 *  \param[in] block_input input matrix to be written in the block
 *  \param[in,out] v_out vector to be updated
 */
template <typename Derived>
inline void setBlock(typename Derived::Index row_min,
                     const Eigen::MatrixBase<Derived> &block_input,
                     Eigen::Ref<SiconosVector> v_out) {
  // Assert if row_min is in range
  assert(row_min < v_out.rows() && "row is out of range");

  // Get the number of rows in the input block
  auto block_rows = block_input.rows();

  // Check if the block fits inside v_out
  assert(row_min + block_rows <= v_out.rows() && "block goes out of range in rows");

  // Using segment to efficiently set the block values in the vector
  v_out.segment(row_min, block_rows) = block_input;
}

/** Set a sub-block of a sparse matrix
 *  \param row_min starting row index of the written block
 *  \param col_min starting col index of the written block
 *  \param[in] block_input input matrix to be written in the block
 *  \param[in,out] m_out matrix to be updated
 */
template <typename Derived>
inline void setBlock(typename Derived::Index row_min, typename Derived::Index col_min,
                     const Eigen::MatrixBase<Derived> &block_input,
                     SiconosSparseMatrix &m_out) {
  // Assert if row and column are in range
  assert(row_min < m_out.rows() && "row is out of range");
  assert(col_min < m_out.cols() && "column is out of range");

  // Get the number of rows and columns in the input block
  auto block_rows = block_input.rows();
  auto block_cols = block_input.cols();

  // Check if the block fits inside m_out
  assert(row_min + block_rows <= m_out.rows() && "block goes out of range in rows");
  assert(col_min + block_cols <= m_out.cols() && "block goes out of range in columns");

  // We need to populate the sparse matrix using triplets
  //  std::vector<Triplet> triplets;
  // Loop through each element in the block_input and add to triplets if non-zero
  for (int i = 0; i < block_rows; ++i) {
    for (int j = 0; j < block_cols; ++j) {
      auto value = block_input(i, j);
      if (value != 0) {  // we can use a tol ?
        m_out.insert(row_min + i, col_min + j) = value;
        // To reset the matrix use
        // triplets.emplace_back(row_min + i, col_min + j, block_input(i, j));
        // + setFromTriplets below
      }
    }
  }
  //  m_out.setFromTriplets(triplets.begin(), triplets.end());
}

}  // namespace siconos::algebra
#endif
