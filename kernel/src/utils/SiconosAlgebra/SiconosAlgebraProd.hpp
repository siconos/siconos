/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#ifndef SA_PROD_HPP
#define SA_PROD_HPP

#include "SiconosAlgebraTypeDef.hpp"
class SiconosMatrix;
class SiconosVector;
class BlockVector;

// NO NEED TO BE FRIEND
/** prod(A, x, y, init) computes y = A*x or y += A*x if init = false
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init = true);

void prod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y, bool init = true);

// NO NEED TO BE FRIEND
/** prod(A, x, y, init) computes y = A*x or y += A*x if init = false
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void prod(const SiconosMatrix& A, const SiconosVector& x, BlockVector& y, bool init = true);

/** prod(x, A, y) computes y = trans(A)*x (init = true) or y += trans(A)*x (init = false)
  \param x a SiconosVector
  \param A a SiconosMatrix
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void prod(const SiconosVector& x, const SiconosMatrix& A, SiconosVector& y, bool init = true);

void prod(const SiconosVector& x, const SiconosMatrix& A, BlockVector& y, bool init = true);

/** prod(A, x) returns the product Ax
  \param A a SiconosMatrix
  \param x a SiconosVector
  \return a SiconosVector
  */
SiconosVector prod(const SiconosMatrix& A, const SiconosVector& x);

/** prod(A, B, C) computes C = A*B in an optimal way, or C += AB if init = false.
  \param A a SiconosMatrix
  \param B a SiconosMatrix
  \param[in,out] C a SiconosMatrix
  \param init a bool (default = true)
  */
void prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C, bool init = true);

/** axpy_prod(A, B, C) computes C = A*B in an optimal way (if init = true, else +=).
  \param A a SiconosMatrix
  \param B a SiconosMatrix
  \param[in,out] C a SiconosMatrix
  \param init a bool (default = true)
  */
void axpy_prod(const SiconosMatrix&, const SiconosMatrix&, SiconosMatrix&, bool);

/** prod(a, A, x, y, init) computes y = a*A*x or y += a*A*x if init = false
  \param a a double
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void prod(double a, const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init = true);

/** subprod(A, x, y) computes sub_y = sub_A*sub_x or sub_y += sub_A*sub_x if init = false
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param coord an Index = [r0A r1A c0A c1A r0x r1x r0y r1y];
  subA is the sub-matrix of A, for row numbers between r0A and r1A-1 and columns between c0A and c1A-1;
  The same for x and y with rix and riy.
  \param init a bool (default = true)
  */
void subprod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, const Index& coord, bool init = true);

void subprod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y, const Index& coord, bool init = true);

/** computes y += sub(transpose(A)) x (only = if init = true) where
subA is a sub-matrix of A.
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param startRow row index of the first element considered in tA (to get sub(tA))
  \param startCol col index of the first element considered in tA (to get sub(tA))
  \param init if true, start with y = 0, else add subA.x to current y.
*/
void taxpy(SPC::SiconosVector x, SPC::SiconosMatrix A, unsigned int startRow, unsigned int startCol, SP::SiconosVector y, bool init = true);

#endif
