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

#include "BlockVector.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace siconos::algebra {

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

void matrixBlockVector_prod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y,
                            bool init = true);

void matrixBlockVector_prod(const SiconosMatrix& A, const BlockVector& x,
                            Eigen::Ref<SiconosVector> y, bool init = true);

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

   \param x input vector
   \param A input matrix
   \param[in,out] y result
   \param init  false to accumulate result into y
  */
void transposeMatrixVector_prod_toBlock(const SiconosVector& x, const SiconosMatrix& A,
                                        BlockVector& y, bool init = true);

}  // namespace siconos::algebra

#endif
