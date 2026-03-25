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

// #include <concepts>
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
bool isSymmetric(const SiconosDenseMatrix& mat, double tol = 1e-12);
bool isSymmetric(const SiconosSparseMatrix& mat, double tol = 1e-12);

// Concept : any eigen matrix complient with middleCols * vec, which means for us dense, sparse
// and map(dense)
template <typename MatrixType>
concept MatrixBlockCompatible =
    requires(const MatrixType& A, siconos::algebra::Index startCol,
             siconos::algebra::Index blockSize, const SiconosVector& vec) {
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
  siconos::algebra::Index startCol = 0;
  for (const auto& it : x) {
    siconos::algebra::Index blockSize = it->size();
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

  siconos::algebra::Index startCol = 0;
  for (auto& block : y) {
    siconos::algebra::Index blockSize = block->size();
    auto subA = A.middleCols(startCol, blockSize);
    *block += subA.transpose() * x;
    startCol += blockSize;
  }
}

/** Generate a random sparse matrix. Useful for tests */
SiconosSparseMatrix generateRandomSparseMatrix(siconos::algebra::Index rows,
                                               siconos::algebra::Index cols,
                                               siconos::algebra::Index nnz,
                                               std::optional<double> density = std::nullopt);

// Tools to solve linear systems (wrap around eigen functions)
// Objective: be able to activate/deactivate residu verification after
// LU solve.

/**  Default signature for 'solving' functions */
using SolveLuFunction = std::function<void(
    const SiconosDenseLUMatrix&, const Eigen::Ref<const SiconosDenseMatrix>&,
    const Eigen::Ref<const SiconosDenseMatrix>&, Eigen::Ref<SiconosDenseMatrix>)>;

/** Default tolerance LU solvers */
inline constexpr double solver_tol = 1e-12;

inline void solve_lu(const SiconosDenseLUMatrix& luw,
                     const Eigen::Ref<const SiconosDenseMatrix>& original,
                     const Eigen::Ref<const SiconosDenseMatrix>& rhs,
                     Eigen::Ref<SiconosDenseMatrix> result) {
  if (result.data() != rhs.data())
    result.noalias() = luw.solve(rhs);
  else
    result = luw.solve(rhs);
}

inline void solve_and_check(const SiconosDenseLUMatrix& luw,
                            const Eigen::Ref<const SiconosDenseMatrix>& original,
                            const Eigen::Ref<const SiconosDenseMatrix>& rhs,
                            Eigen::Ref<SiconosDenseMatrix> result) {
  double residu = 0.;
  if (result.data() != rhs.data()) {
    result.noalias() = luw.solve(rhs);
    // To avoid division by 0
    auto nb = std::max(rhs.norm(), std::sqrt(std::numeric_limits<double>::epsilon()));
    residu = (original * result - rhs).norm() / nb;
  } else {  // we need a copy of the rhs ...
    siconos::algebra::SiconosVector rhs_2 = rhs;
    result.noalias() = luw.solve(rhs_2);

    // if original is not provided. Test purpose.
    //  auto residu = (luw.reconstructedMatrix() * result - rhs).norm();
    // To avoid division by 0
    auto nb = std::max(rhs_2.norm(), std::sqrt(std::numeric_limits<double>::epsilon()));
    residu = (original * result - rhs_2).norm() / nb;
  }
  if (residu > solver_tol) {
    std::cerr << "WARNING: LU solve error exceed the given tolerance. Error=" << residu
              << "\n";
  }
}
/**  Default signature for 'solving' functions */
using SolveInverseFunction =
    std::function<void(const SiconosDenseLUMatrix&,
                       const Eigen::Ref<const SiconosDenseMatrix>&, SiconosDenseMatrix&)>;

inline void inverse_with_check(const SiconosDenseLUMatrix& luw,
                               const Eigen::Ref<const SiconosDenseMatrix>& original,
                               SiconosDenseMatrix& result) {
  result.noalias() = luw.inverse();  // noalias : result != LUw (types différents)

  // Check
  Index n = original.rows();
  double residu = (original * result - SiconosDenseMatrix::Identity(n, n)).norm() / n;
  if (residu > solver_tol)
    std::cerr << "WARNING:: inverse computation results in a residu which exceeds the "
                 "tolerance, error="
              << residu << "\n";
}
inline void inverse_std(const SiconosDenseLUMatrix& luw,
                        const Eigen::Ref<const SiconosDenseMatrix>& original,
                        SiconosDenseMatrix& result) {
  result.noalias() = luw.inverse();  // noalias : result != LUw (types différents)
}

/** @brief solve Ax=b with LU decomposition
 *  @param luw LU factorisation of the matrix
 *  @param original A matrix, not factorized. Unused with "no check" mode
 *  @param rhs righ-hand side
 *  @param result solution
 */
inline SolveLuFunction solve = solve_lu;

/** @brief compute the inverse of a matrix (from LU dec)
 *  @param luw LU factorisation of the matrix
 *  @param original matrix, not factorized. Unused with "no check" mode
 *  @param result inverse
 */
inline SolveInverseFunction inverse = inverse_std;

/** True to activate result check in linear system solve */
inline bool is_solver_check_activated = false;

/** Activate solver checking */
inline void enable_solver_check() {
  is_solver_check_activated = true;
  solve = solve_and_check;
  inverse = inverse_with_check;
}

/** De-activate solver checking */
inline void disable_solver_check() {
  is_solver_check_activated = false;
  solve = solve_lu;
  inverse = inverse_std;
}

}  // namespace siconos::algebra

#endif
