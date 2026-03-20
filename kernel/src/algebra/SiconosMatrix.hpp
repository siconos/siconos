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

// #define EIGEN_DEFAULT_DENSE_INDEX_TYPE long long // Done using CMake def.

#include "EigenInclude.hpp"  // IWYU pragma: keep - Must be included before Eigen/Core
//
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/SparseCore>  // For Eigen Sparse matrices
#include <Eigen/SparseLU>

#include "CSparseMatrix.h"  // IWYU pragma: keep
#include "SiconosVector.hpp"

// Default type for sparse storage indices
// Rq: dense storage indices type is defined by
// EIGEN_DEFAULT_DENSE_INDEX_TYPE
// Both SPARSE_STORAGE_INDEX and EIGEN_DEFAULT_DENSE_INDEX_TYPE
// are supposed to be set during cmake conf
// (look for target_compile_definition ...)
#ifndef SICONOS_SPARSE_STORAGE_INDEX_TYPE
#define SICONOS_SPARSE_STORAGE_INDEX_TYPE EIGEN_DEFAULT_DENSE_INDEX_TYPE
#endif

struct NumericsMatrix;

namespace siconos::algebra {

/**
   Abstract class to provide interface for matrices handling

   Matrices can be either block or Simple.
   See Derived classes for details.

   You can find an overview on how to build and use vectors and matrices in
   siconos users guide
   .

*/
using SiconosDenseMatrix =
    Eigen::Matrix<double_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using SiconosMatrix33 = Eigen::Matrix<double_t, 3, 3, Eigen::ColMajor>;
using SiconosMatrix36 = Eigen::Matrix<double_t, 3, 6, Eigen::ColMajor>;
using SiconosMatrix37 = Eigen::Matrix<double_t, 3, 7, Eigen::ColMajor>;
using SiconosMatrix66 = Eigen::Matrix<double_t, 6, 6, Eigen::ColMajor>;
using SiconosMatrix76 = Eigen::Matrix<double_t, 7, 6, Eigen::ColMajor>;
using SiconosMatrix67 = Eigen::Matrix<double_t, 6, 7, Eigen::ColMajor>;
using SiconosMatrix33Diagonal = Eigen::DiagonalMatrix<double_t, 3>;
using SiconosMatrix66Diagonal = Eigen::DiagonalMatrix<double_t, 6>;

inline const SiconosMatrix33 identity33 = SiconosMatrix33::Identity();

/** Sparse matrix storage */
using SparseIndex = SICONOS_SPARSE_STORAGE_INDEX_TYPE;
using SiconosSparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, SparseIndex>;
using Triplet = Eigen::Triplet<double>;  // Used to fill sparse matrices

/** Utility function to ensure proper cast from Sparse Eigen indices type to size_t (or
 * equivalent) */
template <std::unsigned_integral T>  // Nodiscard to force a return value (else warning)
[[nodiscard]] T sparse_index_to_unsigned(SparseIndex value) {
  // if (!std::in_range<T>(value))
  //   throw std::overflow_error(
  //       "SparseIndex value cannot be represented as an unsigned integer");
  assert(std::in_range<T>(value));
  return static_cast<T>(value);
}

// Map types

using SiconosDiagonalMatrix = Eigen::DiagonalMatrix<double_t, Eigen::Dynamic>;

// Tmp version. To replace BlockVector ?
using EigenBlock = std::vector<Eigen::Ref<SiconosVector>>;

// LU Factorization
using SiconosDenseLUMatrix = Eigen::FullPivLU<siconos::algebra::SiconosDenseMatrix>;
using SiconosSparseLUMatrix = Eigen::SparseLU<SiconosSparseMatrix>;

// --- Select between dense and sparse storage ---
// SICONOS_SPARSE is set by cmake.
// Use -DSICONOS_SPARSE=1 to activate sparse storage
#ifdef SICONOS_SPARSE
using SiconosMatrix = SiconosSparseMatrix;
using SiconosLUMatrix = SiconosSparseLUMatrix;
#else
using SiconosMatrix = SiconosDenseMatrix;
using SiconosLUMatrix = SiconosDenseLUMatrix;
#endif

using MapType = Eigen::Map<SiconosDenseMatrix>;
using ConstMapType = Eigen::Map<const SiconosDenseMatrix>;

enum class StorageType { dense, sparse };

/** General dense matrix print function
 *  \param mat the matrix to be displayed
 */
template <typename Derived>
void print(const Eigen::EigenBase<Derived>& mat) {
  std::cout << std::scientific << std::setprecision(6) << mat.derived() << std::endl;
}

// Specialization for sparse matrices
inline void print(const SiconosSparseMatrix& mat) {
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
// inline void print(const SiconosVector& vec) {
//   std::cout << std::scientific << std::setprecision(6) << vec.transpose() << "\n";
// }
template <typename Derived>
  requires(Derived::ColsAtCompileTime == 1 || Derived::RowsAtCompileTime == 1)
void print(const Eigen::MatrixBase<Derived>& vec) {
  std::cout << std::scientific << std::setprecision(6) << vec.transpose() << "\n";
}

/** \return inf. norm of the matrix
 *  \param mat the input matrix
 */
template <typename Derived>
auto normInf(const Eigen::MatrixBase<Derived>& mat) {
  return mat.cwiseAbs().rowwise().sum().maxCoeff();
}

// Specialization for sparse matrices (uses iterators)
inline double normInf(const SiconosSparseMatrix& mat) {
  // Warn: col-major
  siconos::algebra::SiconosVector rowSums{mat.rows()};
  rowSums.setZero();
  for (int k = 0; k < mat.outerSize(); ++k) {
    // Sum over each row
    for (SiconosSparseMatrix::InnerIterator it(mat, k); it; ++it) {
      rowSums(it.row()) += std::abs(it.value());
    }
  }
  return rowSums.maxCoeff();
}

/** Set a sub-block of a matrix
 *  \param row_min starting row index of the written block
 *  \param col_min starting col index of the written block
 *  \param[in] input input matrix to be written in the block
 *  \param[in,out] m_out matrix to be updated
 */
template <typename Derived>
inline void setBlock(typename Derived::Index row_min, typename Derived::Index col_min,
                     const Eigen::MatrixBase<Derived>& input,
                     Eigen::Ref<SiconosDenseMatrix> m_out) {
  // assert(m != *this);
  assert(row_min < m_out.rows() && "row is out of range");
  assert(col_min < m_out.cols() && "column is out of range");

  // Get the number of rows and columns in the input block
  auto block_rows = input.rows();
  auto block_cols = input.cols();

  // Calculate the max row and column indices
  assert(row_min + block_rows <= m_out.rows() && "block goes out of range in rows");
  assert(col_min + block_cols <= m_out.cols() && "block goes out of range in columns");

  m_out.block(row_min, col_min, block_rows, block_cols) = input;
}

// /** Set a sub-block of a Vector
//  *  \param row_min starting row index of the written block
//  *  \param[in] block_input input matrix to be written in the block
//  *  \param[in,out] v_out vector to be updated
//  */
// template <typename Derived>
// inline void setBlock(typename Derived::Index row_min,
//                      const Eigen::MatrixBase<Derived> &block_input,
//                      Eigen::Ref<SiconosVector> v_out) {
//   // Assert if row_min is in range
//   assert(row_min < v_out.rows() && "row is out of range");

//   // Get the number of rows in the input block
//   auto block_rows = block_input.rows();

//   // Check if the block fits inside v_out
//   assert(row_min + block_rows <= v_out.rows() && "block goes out of range in rows");

//   // Using segment to efficiently set the block values in the vector
//   v_out.segment(row_min, block_rows) = block_input;
// }

/** Set a sub-block of a sparse matrix
 *  \param row_min starting row index of the written block
 *  \param col_min starting col index of the written block
 *  \param[in] input input matrix to be written in the block
 *  \param[in,out] m_out matrix to be updated
 */
template <typename Derived>
inline void setBlock(typename Derived::Index row_min, typename Derived::Index col_min,
                     const Eigen::MatrixBase<Derived>& input, SiconosSparseMatrix& m_out) {
  // Assert if row and column are in range
  assert(row_min < m_out.rows() && "row is out of range");
  assert(col_min < m_out.cols() && "column is out of range");

  // Get the number of rows and columns in the input block
  auto block_rows = input.rows();
  auto block_cols = input.cols();

  // Check if the block fits inside m_out
  assert(row_min + block_rows <= m_out.rows() && "block goes out of range in rows");
  assert(col_min + block_cols <= m_out.cols() && "block goes out of range in columns");

  // We need to populate the sparse matrix using triplets
  //  std::vector<Triplet> triplets;
  // Loop through each element in the input and add to triplets if non-zero
  if (m_out.isCompressed()) m_out.uncompress();
  for (int i = 0; i < block_rows; ++i) {
    for (int j = 0; j < block_cols; ++j) {
      auto value = input(i, j);
      if (value != 0) {  // we can use a tol ?
        m_out.coeffRef(row_min + i, col_min + j) = value;
        // To reset the matrix use
        // triplets.emplace_back(row_min + i, col_min + j, input(i, j));
        // + setFromTriplets below
      }
    }
  }
  m_out.makeCompressed();
  //  m_out.setFromTriplets(triplets.begin(), triplets.end());
}

}  // namespace siconos::algebra
#endif
