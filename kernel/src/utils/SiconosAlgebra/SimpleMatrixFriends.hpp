/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/*! \file SimpleMatrixFriends.hpp
  List of friend functions for SimpleMatrix.
  */

#ifndef __SimpleMatrixFriends__
#define __SimpleMatrixFriends__
#include "SiconosMatrix.hpp"

/** Copy a subBlock of MIn into a sub-block of MOut - Dim and positions of the sub-block are given in dim and start.
 * \param MIn a SPC::SiconosMatrix
 * \param[in,out] MOut a SP::SiconosMatrix
 * \param dim an Index, dim[0], dim[1]: number of rows and columns of the sub-block
 * \param start an Index, start[0], start[1]: position (row, column) of the first element of the sub-block in MIn
 *  start[2], start[3]: position (row, column) of the first element of the sub-block in MOut.
 */
void setBlock(SPC::SiconosMatrix MIn, SP::SiconosMatrix MOut, const Index& dim, const Index& start);

/** multiplication of a matrix by a double
 *  \param A a SiconosMatrix
 *  \param a a double
 *  \return a SimpleMatrix
 */
const SimpleMatrix operator * (const SiconosMatrix& A, double a);

/** multiplication of a SimpleMatrix by a SimpleMatrix
 *  \param a const SiconosMatrix&
 *  \param a const SimpleMatrix&
 *  \return a const SimpleMatrix
 */
//  const SimpleMatrix operator * (const SimpleMatrix&,const SimpleMatrix&);

/** multiplication of a SP::SimpleMatrix by a SP::SimpleMatrix
 *  \param A a SP::SiconosMatrix
 *  \param B a SP::SimpleMatrix
 *  \return a SP::SimpleMatrix
 */
SP::SimpleMatrix operator * (const SP::SimpleMatrix A, const SP::SimpleMatrix B);

/** operator += add B to A
 *  \param[in,out] A a SP::SiconosMatrix
 *  \param B a SP::SiconosMatrix
 */
void operator +=(SP::SiconosMatrix A, SP::SimpleMatrix B);

/** multiplication of a matrix by a double
 *  \param a a double
 *  \param A a SiconosMatrix
 *  \return a SimpleMatrix
 */
SimpleMatrix operator * (double a, const SiconosMatrix& A);

/** division of the matrix by a double
 *  \param A a SiconosMatrix
 *  \param a a double
 *  \return a SimpleMatrix
 */
const SimpleMatrix operator /(const SiconosMatrix& A, double a);

/** Addition of two matrices, C = A+B
 * \param A a SiconosMatrix
 * \param B a SiconosMatrix
 * \return a SimpleMatrix C
 */
const SimpleMatrix operator +(const SiconosMatrix& A, const SiconosMatrix& B);
/** Addition of two matrices, C = A+B
 * \param A a SP::SiconosMatrix
 * \param B a SP::SiconosMatrix
 * \return a SP::SimpleMatrix
 */
SP::SimpleMatrix operator +(const SP::SimpleMatrix A, const SP::SimpleMatrix B);
SimpleMatrix operator +(const SimpleMatrix& A, const SimpleMatrix& B);

/** Addition of two matrices C = A+B
 *  \param A a SiconosMatrix
 *  \param B a SiconosMatrix
 *  \param[in,out] C a SiconosMatrix
 */
void add(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C);

/** Subtraction of two matrices, C = A-B
 * \param A a SiconosMatrix
 * \param B a SiconosMatrix
 * \return a SimpleMatrix
 */
const SimpleMatrix operator -(const SiconosMatrix& A, const SiconosMatrix& B);
//  const SimpleMatrix operator -(const SimpleMatrix&,const SimpleMatrix&);

/** Subtraction of two matrices C = A-B
 *  \param A a SiconosMatrix
 *  \param B a SiconosMatrix
 *  \param[in,out] C a SiconosMatrix
 */
void sub(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C);

/**: A==B when (A-B).normInf()<tolerance
   \param A a SiconosMatrix
   \param B a SiconosMatrix
   \return a boolean
 */
bool operator == (const SiconosMatrix& A, const SiconosMatrix& B);

/**: A!=B when (A-B).normInf()>tolerance
 * \param A a SiconosMatrix
 * \param B a SiconosMatrix
 * \return a boolean
 */
bool operator != (const SiconosMatrix& A, const SiconosMatrix& B);

/** compute the power of the matrix (!)
 * \param A a SimpleMatrix
 * \param e the exponent (an unsigned int)
 * \return a SimpleMatrix
 */
const SimpleMatrix matrix_pow(const SimpleMatrix& A, unsigned int e);

/** product of two matrices, C = A*B
  \param A a SiconosMatrix
  \param B a SiconosMatrix
  \return C a SimpleMatrix
  */
const SimpleMatrix prod(const SiconosMatrix& A, const SiconosMatrix& B);

/** prod(A, B, C) computes C = A*B in an optimal way, or C += AB if init = false.
  \param A a SiconosMatrix
  \param B a SiconosMatrix
  \param[in,out] C a SiconosMatrix
  \param init a bool (default = true)
  */
void prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C, bool init = true);

/** prod(A, B, C) computes C = A*B in an optimal way (if init = true, else +=).
  \param A a SiconosMatrix
  \param B a SiconosMatrix
  \param[in,out] C a SiconosMatrix
  \param init a bool (default = true)
  */
void axpy_prod(const SiconosMatrix&, const SiconosMatrix&, SiconosMatrix&, bool);

/** prod(A, x) returns the product Ax
  \param A a SiconosMatrix
  \param x a SiconosVector
  \return a SiconosVector
  */
const SiconosVector prod(const SiconosMatrix& A, const SiconosVector& x);

/** prod(A, x, y, init) computes y = A*x or y += A*x if init = false
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init = true);

/** prod(A, x, y, init) computes y = A*x or y += A*x if init = false
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void prod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y, bool init = true);

/** prod(A, x, y, init) computes y = A*x or y += A*x if init = false
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void prod(const SiconosMatrix& A, const SiconosVector& x, BlockVector& y, bool init = true);

/** prod(a, A, x, y, init) computes y = a*A*x or y += a*A*x if init = false
  \param a a double
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void prod(double a, const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init = true);

/** prod(x, A, y) computes y = trans(A)*x (init = true) or y += trans(A)*x (init = false)
  \param x a SiconosVector
  \param A a SiconosMatrix
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void prod(const SiconosVector& x, const SiconosMatrix& A, SiconosVector& y, bool init = true);
void prod(const SiconosVector& x, const SiconosMatrix& A, BlockVector& y, bool init = true);

/** subprod(A, x, y) computes sub_y = sub_A*sub_x or sub_y += sub_A*sub_x if init = false
  \param A a SiconosMatrix
  \param x a BlockVector
  \param[in,out] y a SiconosVector
  \param coord an Index = [r0A r1A c0A c1A r0x r1x r0y r1y];
  subA is the sub-matrix of A, for row numbers between r0A and r1A-1 and columns between c0A and c1A-1;
  The same for x and y with rix and riy.
  \param init a bool (default = true)
  */
void subprod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y, const Index& coord, bool init = true);

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

/** computes y = A*x (init = true) or y += A*x (init = false)
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param init a bool
  */
void axpy_prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init);

/** gemvtranspose(a, A, x, b, y) computes y = a*transpose(A)*x + b*y
    This function wrapped blas gemv through boost bindings. Reserved to dense matrices and vectors.
    \param a a double
    \param A a SiconosMatrix
    \param x a SiconosVector
    \param b a double
    \param[in,out] y a SiconosVector
*/
void gemvtranspose(double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y);

/** gemv(a, A, x, b, y) computes y = a*A*x+ b*y
    This function wrapped blas gemv through boost bindings. Reserved to dense matrices and vectors.
    \param a a double
    \param A a SiconosMatrix
    \param x a SiconosVector
    \param b a double
    \param[in,out] y a SiconosVector
*/
void gemv(double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y);

/** gemmtranspose(a, A, B, b, C) computes C = a*transpose(A)*transpose(B) + b*C
    This function wrapped blas gemm through boost bindings. Reserved to dense matrices and vectors.
    \param a a double
    \param A a SiconosMatrix
    \param B a SiconosMatrix
    \param b a double
    \param[in,out] C a SiconosMatrix
*/
void gemmtranspose(double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C);

/** gemm(a, A, B, b, C) computes C = a*A*B+ b*C
    This function wrapped blas gemm through boost bindings. Reserved to dense matrices and vectors.
    \param a a double
    \param A a SiconosMatrix
    \param B a SiconosMatrix
    \param b a double
    \param[in,out] C a SiconosMatrix
*/
void gemm(double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C);

/** multiplication of a matrix by a scalar, B = a*A (init = true) or B += a*A (init = false)
 *  \param a a double
 *  \param A a SiconosMatrix
 *  \param[in,out] B a SiconosMatrix
 *  \param init a bool
 */
void scal(double a, const SiconosMatrix& A, SiconosMatrix& B, bool = true);

void private_addprod(const SiconosMatrix& , unsigned int, unsigned int, const SiconosVector&, SiconosVector&);
void private_addprod(double, SPC::SiconosMatrix, unsigned int, unsigned int, SPC::SiconosVector, SP::SiconosVector);
void private_addprod(SPC::SiconosVector, SPC::SiconosMatrix, unsigned int, unsigned int, SP::SiconosVector);
void private_addprod(const SiconosMatrix& , unsigned int, unsigned int, const BlockVector&, SiconosVector&);
void private_addprod(SPC::BlockVector, SPC::SiconosMatrix, unsigned int, unsigned int, SP::SiconosVector);
void private_prod(const SiconosMatrix& A, unsigned int, const SiconosVector&, SiconosVector&, bool);
void private_prod(const SiconosMatrix& A, unsigned int, const BlockVector& , SiconosVector&, bool);
void private_prod(SPC::SiconosMatrix, unsigned int, SPC::SiconosVector , SP::BlockVector, bool);
void private_prod(SPC::SiconosMatrix, unsigned int, SPC::BlockVector , SP::BlockVector, bool);
void private_prod(double, SPC::SiconosMatrix, unsigned int, SPC::SiconosVector, SP::SiconosVector, bool);
void private_prod(SPC::SiconosVector, SPC::SiconosMatrix, unsigned int, SP::SiconosVector, bool);
void private_prod(SPC::BlockVector, SPC::SiconosMatrix, unsigned int, SP::SiconosVector, bool);
void private_prod(SPC::BlockVector, SPC::SiconosMatrix, unsigned int, SP::BlockVector, bool);
void private_prod(SPC::SiconosVector, SPC::SiconosMatrix, unsigned int, SP::BlockVector, bool);

#endif
