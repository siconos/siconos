/* Siconos-Kernel, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*! \file SimpleMatrixFriends.hpp
  List of friend functions for SimpleMatrix.
  */

#ifndef __SimpleMatrixFriends__
#define __SimpleMatrixFriends__
#include "SiconosMatrix.hpp"
//#include "BlocksSimpleMat.hpp"

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

/**
 *Default comparator
 *\param A a SimpleMatrix
 *\param B a SimpleMatrix
 * return true if A != B
 */
bool operator!= (const SimpleMatrix& A, const SimpleMatrix& B);

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
//  SimpleMatrix operator +(const SimpleMatrix&,const SimpleMatrix&);

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
 * \param A a SiconosMatrix
 * \param B a SiconosMatrix
 * \return a boolean
 */
bool operator == (const SiconosMatrix& A, const SiconosMatrix& B);

/** compute the power of the matrix (!)
 * \param A a SimpleMatrix
 * \param e the exponent (an unsigned int)
 * \return a SimpleMatrix
 */
const SimpleMatrix pow(const SimpleMatrix& A, unsigned int e);

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

/** prod(A, B, C) computes C = A*B in an optimal way (= if init = true, else +=).
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

/** prod(A, x, y) computes y = A*x, using atlas::gemv - Reserved to dense matrices and vectors.
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  */
void gemv(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y);

/** gemv(transA, a, A, x, b, y) computes y = a*op(A)*x + b*y
  with transX = CblasNoTrans (op(X) = X), CblasTrans (op(X) = transpose(X)), CblasConjTrans (op(X) = conj(X))
  This function encapsulates atlas::gemv. Reserved to dense matrices and vectors.
  \param transA a CBLAS_TRANSPOSE, op for A
  \param a a double
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param b a double
  \param[in,out] y a SiconosVector
  */
void gemv(const CBLAS_TRANSPOSE transA, double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y);

/** gemv(a, A, x, b, y) computes y = a*A*x+ b*y
  This function encapsulates atlas::gemv. Reserved to dense matrices and vectors.
  \param a a double
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param b a double
  \param[in,out] y a SiconosVector
  */
void gemv(double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y);

/** gemm(transA, transB, a, A, B, b, C) computes C = a*op(A)*op(B) + b*C
  with transX = CblasNoTrans (op(X) = X), CblasTrans (op(X) = transpose(X)), CblasConjTrans (op(X) = conj(X))
  This function encapsulates atlas::gemm. Reserved to dense matrices.
  \param transA a CBLAS_TRANSPOSE, op for A
  \param transB a CBLAS_TRANSPOSE, op for B
  \param a a double
  \param A a SiconosMatrix
  \param B a SiconosMatrix
  \param b a double
  \param[in,out] C a SiconosMatrix
  */
void gemm(const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB, double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C);

/** gemm(a, A, B, b, C) computes C = a*A*B+ b*C
  This function encapsulates atlas::gemm. Reserved to dense matrices.
  \param a a double
  \param A a SiconosMatrix
  \param B a SiconosMatrix
  \param b a double
  \param[in,out] C a SiconosMatrix
  */
void gemm(double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C);

/** gemm(A, B, C) computes C = A*B
  This function encapsulates atlas::gemm. Reserved to dense matrices.
  \param A a SiconosMatrix
  \param B a SiconosMatrix
  \param[in,out] C a SiconosMatrix
  */
void gemm(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C);

/** multiplication of a matrix by a scalar, B = a*A (init = true) or B += a*A (init = false)
 *  \param a a double
 *  \param A a SiconosMatrix
 *  \param[in,out] B a SiconosMatrix
 *  \param init a bool
 */
void scal(double a, const SiconosMatrix& A, SiconosMatrix& B, bool = true);

#endif
