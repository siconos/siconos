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
using SiconosMatrix = Eigen::Matrix<double_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using SiconosMatrix33 = Eigen::Matrix<double_t, 3, 3, Eigen::ColMajor>;
/** Sparse matrix storage */
using SiconosSparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
using Triplet = Eigen::Triplet<double>;  // Used to fill sparse matrices

// Map types

using MapVectorType = Eigen::Map<SiconosVector>;
using MapVector3Type = Eigen::Map<SiconosVector3>;
using ConstMapVectorType = Eigen::Map<const SiconosVector>;
using ConstMapVector3Type = Eigen::Map<const SiconosVector3>;
using MapType = Eigen::Map<SiconosMatrix>;
using ConstMapType = Eigen::Map<const SiconosMatrix>;

using SiconosDiagonalMatrix = Eigen::DiagonalMatrix<double_t, Eigen::Dynamic>;

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

}  // namespace siconos::algebra
#endif
