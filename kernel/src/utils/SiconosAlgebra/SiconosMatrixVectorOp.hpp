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

/*! \file SiconosMatrixVectorOp.hpp
  \brief Toolbox for operators acting on matrices and vectors
*/

#ifndef SICOMATVEC_OPH
#define SICOMATVEC_OPH

// #include <boost/numeric/ublas/fwd.hpp>
#include <complex>

#include "BlockVector.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace siconos::algebra {

/** prod(A, x, y, init) computes y = A*x or y += A*x if init = false
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init = true);
// void prod(const MapType& A, const SiconosVector& x, SiconosVector& y, bool init = true);

void matrixBlockVector_prod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y,
                            bool init = true);

// NO NEED TO BE FRIEND
/** computes y = A*x or y += A*x if init = false
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void matrixVector_prod_toBlock(const SiconosMatrix& A, const SiconosVector& x, BlockVector& y,
                               bool init = true);

/** computes y = trans(A)*x (init = true) or y += trans(A)*x (init = false)
  \param x a SiconosVector
  \param A a SiconosMatrix
  \param[in,out] y a SiconosVector
  \param init false to accumulate result into y
  */
void transposeMatrixVector_prod(const SiconosVector& x, const SiconosMatrix& A,
                                SiconosVector& y, bool init = true);

/** computes y(block vector) = trans(A)*x (init = true) or y += trans(A)*x (init= false)

   \param x input vector
   \param A input matrix
   \param[in,out] y result
   \param init  false to accumulate result into y
  */
void transposeMatrixVector_prod_toBlock(const SiconosVector& x, const SiconosMatrix& A,
                                        BlockVector& y, bool init = true);

/** prod(a, A, x, y, init) computes y = a*A*x or y += a*A*x if init = false
  \param a a double
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void prod(double a, const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y,
          bool init = true);

/** computes y += sub(transpose(A)) x (only = if init = true) where
subA is a sub-matrix of A.
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param startRow row index of the first element considered in tA (to get sub(tA))
  \param startCol col index of the first element considered in tA (to get sub(tA))
  \param init if true, start with y = 0, else add subA.x to current y.
*/
void taxpy(const SiconosVector& x, const SiconosMatrix& A, unsigned int startRow,
           unsigned int startCol, std::shared_ptr<SiconosVector> y, bool init = true);

}  // namespace siconos::algebra

#endif
