/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

DEFINE_SPTR(SimpleMatrix);

/** Copy a subBlock of MIn into a sub-block of MOut - Dim and positions of the sub-block are given in dim and start.
 * \param MIn, a SiconosMatrix*
 * \param MOut, a SiconosMatrix*
 * \param dim, an Index, dim[0], dim[1]: number of rows and columns of the sub-block
 * \param start, an Index, start[0], start[1]: position (row, column) of the first element of the sub-block in MIn
 *  start[2], start[3]: position (row, column) of the first element of the sub-block in MOut.
 */
void setBlock(SPC::SiconosMatrix , SP::SiconosMatrix , const Index&, const Index&);
/** multiplication of a matrix by a double
 *  \param a SiconosMatrix
 *  \param a double
 *  \return a SimpleMatrix
 */
const SimpleMatrix operator * (const SiconosMatrix&, double);
/** multiplication of a SimpleMatrix by a SimpleMatrix
 *  \param a const SiconosMatrix&
 *  \param a const SimpleMatrix&
 *  \return a const SimpleMatrix
 */
//  const SimpleMatrix operator * (const SimpleMatrix&,const SimpleMatrix&);
/** multiplication of a SP::SimpleMatrix by a SP::SimpleMatrix
 *  \param a const SP::SiconosMatrix
 *  \param a const SP::SimpleMatrix
 *  \return a const SP::SimpleMatrix
 */
SP::SimpleMatrix operator * (const SP::SimpleMatrix, const SP::SimpleMatrix);

/**
 *Default comparator
 *\param a SP::SiconosMatrix
 *\param a SP::SiconosMatrix
 * return true
 */
bool operator!= (const SimpleMatrix&, const SimpleMatrix&);

/** operator += add SP::SimpleMatrix
 *  \param in/outSP::SiconosMatrix : a matrix A
 *  \param SP::SiconosMatrix : a matrix B
 */
void operator +=(SP::SiconosMatrix, SP::SimpleMatrix);

/** multiplication of a matrix by a double
 *  \param a double
 *  \param a SiconosMatrix
 *  \return a SimpleMatrix
 */
SimpleMatrix operator * (double , const SiconosMatrix&);

/** division of the matrix by a double
 *  \param a SiconosMatrix
 *  \param a double
 *  \return a SimpleMatrix
 */
const SimpleMatrix operator /(const SiconosMatrix&, double);

/** Addition of two matrices, C = A+B
 * \param SiconosMatrix A
 * \param SiconosMatrix B
 * \return a SimpleMatrix C
 */
const SimpleMatrix operator +(const SiconosMatrix&, const SiconosMatrix&);
/** Addition of two matrices, C = A+B
 * \param SP::SiconosMatrix A
 * \param SP::SiconosMatrix B
 * \return a SP::SimpleMatrix
 */
SP::SimpleMatrix operator +(const SP::SimpleMatrix, const SP::SimpleMatrix);
//  SimpleMatrix operator +(const SimpleMatrix&,const SimpleMatrix&);

/** Addition of two matrices C = A+B
 *  \param SiconosMatrix A (in)
 *  \param SiconosMatrix B (in)
 *  \param SiconosMatrix C (in-out)
 */
void add(const SiconosMatrix&, const SiconosMatrix&, SiconosMatrix&);

/** Subtraction of two matrices, C = A-B
 * \param SiconosMatrix A
 * \param SiconosMatrix B
 * \return a SimpleMatrix C
 */
const SimpleMatrix operator -(const SiconosMatrix&, const SiconosMatrix&);
//  const SimpleMatrix operator -(const SimpleMatrix&,const SimpleMatrix&);

/** Subtraction of two matrices C = A-B
 *  \param SiconosMatrix A (in)
 *  \param SiconosMatrix B (in)
 *  \param SiconosMatrix C (in-out)
 */
void sub(const SiconosMatrix&, const SiconosMatrix&, SiconosMatrix&);

/**: A==B when (A-B).normInf()<tolerance
 * \param SiconosMatrix A
 * \param SiconosMatrix B
 * \return a boolean
 */
bool operator == (const SiconosMatrix&, const SiconosMatrix&);

/*   /\** compute the power of the matrix (!) */
/*    *  \return a SimpleMatrix */
/*    *\/ */
const SimpleMatrix pow(const SimpleMatrix&, unsigned int);

/** product of two matrices, C = A*B
 *  \param SiconosMatrix A (in)
 *  \param SiconosMatrix B (in)
 *  \return SimpleMatrix C
 \param init, a bool (default = true)
*/
const SimpleMatrix prod(const SiconosMatrix&, const SiconosMatrix&);

/** prod(A, B, C) computes C = A*B in an optimal way, or C += AB if init = false.
    \param a SiconosMatrix, A (in)
    \param a SiconosMatrix, B (in)
    \param a SiconosMatrix, C (in-out)
    \param init, a bool (default = true)
*/
void prod(const SiconosMatrix&, const SiconosMatrix&, SiconosMatrix&, bool = true);

/** prod(A, B, C) computes C = A*B in an optimal way (= if init = true, else +=).
    \param A, a SiconosMatrix
    \param B, a SiconosMatrix
    \param C, a SiconosMatrix
    \param init, a bool
*/
void axpy_prod(const SiconosMatrix&, const SiconosMatrix&, SiconosMatrix&, bool);

/*   /\** compute the product matrix-vector */
/*    *  \return a SimpleVector */
/*    *\/ */
const SimpleVector prod(const SiconosMatrix&, const SiconosVector&);

/** prod(A, x, y) computes y = A*x or y += A*x if init = false
    \param a SiconosMatrix, A (in)
    \param a SiconosVector, x (in)
    \param a SiconosVector, y (in-out)
    \param init, a bool (default = true)
*/
void prod(const SiconosMatrix&, const SiconosVector&, SiconosVector&, bool = true);

/** prod(a, A, x, y, init) computes y = a*A*x or y += a*A*x if init = false
    \param double, a (in)
    \param a SiconosMatrix, A (in)
    \param a SiconosVector, x (in)
    \param a SiconosVector, y (in-out)
    \param init, a bool (default = true)
*/
void prod(double, const SiconosMatrix&, const SiconosVector&, SiconosVector&, bool = true);

/** prod(x, A, y) computes y = trans(A)*x (init = true) or y += trans(A)*x (init = false)
    \param a SiconosVector, x (in)
    \param a SiconosMatrix, A (in)
    \param a SiconosVector, y (in-out)
    \param a bool, init
*/
void prod(const SiconosVector&, const SiconosMatrix&, SiconosVector&, bool = true);

/** subprod(A, x, y) computes sub_y = sub_A*sub_x or sub_y += sub_A*sub_x if init = false
    \param a SiconosMatrix, A (in)
    \param a SiconosVector, x (in)
    \param a SiconosVector, y (in-out)
    \param an Index = [r0A r1A c0A c1A r0x r1x r0y r1y];
    subA is the sub-matrix of A, for row numbers between r0A and r1A-1 and columns between c0A and c1A-1;
    The same for x and y with rix and riy.
    \param init, a bool (default = true)
*/
void subprod(const SiconosMatrix&, const SiconosVector&, SiconosVector&, const Index&, bool = true);

/** computes y = A*x (init = true) or y += A*x (init = false)
    \param a SiconosMatrix, A (in)
    \param a SiconosVector, x (in)
    \param a SiconosVector, y (in-out)
    \param init, a bool
*/
void axpy_prod(const SiconosMatrix&, const SiconosVector&, SiconosVector&, bool);

/** prod(A, x, y) computes y = A*x, using atlas::gemv - Reserved to dense matrices and vectors.
    \param a SiconosMatrix, A (in)
    \param a SiconosVector, x (in)
    \param a SiconosVector, y (in-out)
*/
void gemv(const SiconosMatrix&, const SiconosVector&, SiconosVector&);

/** gemv(transA, a, A, x, b, y) computes y = a*op(A)*x + b*y
    with transX = CblasNoTrans (op(X) = X), CblasTrans (op(X) = transpose(X)), CblasConjTrans (op(X) = conj(X))
    This function encapsulates atlas::gemv. Reserved to dense matrices and vectors.
    \param CBLAS_TRANSPOSE, op for A
    \param a double, a (in)
    \param a SiconosMatrix, A (in)
    \param a SiconosVector, x (in)
    \param a double, b (in)
    \param a SiconosVector, y (in-out)
*/
void gemv(CBLAS_TRANSPOSE, double, const SiconosMatrix&, const SiconosVector&, double, SiconosVector&);

/** gemv(a, A, x, b, y) computes y = a*A*x+ b*y
    This function encapsulates atlas::gemv. Reserved to dense matrices and vectors.
    \param a double, a (in)
    \param a SiconosMatrix, A (in)
    \param a SiconosVector, x (in)
    \param a double, b (in)
    \param a SiconosVector, y (in-out)
*/
void gemv(double, const SiconosMatrix&, const SiconosVector&, double, SiconosVector&);

/** gemm(transA, transB, a, A, B, b, C) computes C = a*op(A)*op(B) + b*C
    with transX = CblasNoTrans (op(X) = X), CblasTrans (op(X) = transpose(X)), CblasConjTrans (op(X) = conj(X))
    This function encapsulates atlas::gemm. Reserved to dense matrices.
    \param CBLAS_TRANSPOSE, op for A
    \param CBLAS_TRANSPOSE, op for B
    \param a double, a (in)
    \param a SiconosMatrix, A (in)
    \param a SiconosMatrix, B (in)
    \param a double, b (in)
    \param a SiconosMatrix, C (in-out)
*/
void gemm(const CBLAS_TRANSPOSE, const CBLAS_TRANSPOSE, double, const SiconosMatrix&, const SiconosMatrix&, double, SiconosMatrix&);

/** gemm(a, A, B, b, C) computes C = a*A*B+ b*C
    This function encapsulates atlas::gemm. Reserved to dense matrices.
    \param a double, a (in)
    \param a SiconosMatrix, A (in)
    \param a SiconosMatrix, B (in)
    \param a double, b (in)
    \param a SiconosMatrix, C (in-out)
*/
void gemm(double, const SiconosMatrix&, const SiconosMatrix&, double, SiconosMatrix&);

/** gemm(A, B, C) computes C = A*B
    This function encapsulates atlas::gemm. Reserved to dense matrices.
    \param a SiconosMatrix, A (in)
    \param a SiconosMatrix, B (in)
    \param a SiconosMatrix, C (in-out)
*/
void gemm(const SiconosMatrix&, const SiconosMatrix&, SiconosMatrix&);

/** multiplication of a matrix by a scalar, B = a*A (init = true) or B += a*A (init = false)
 *  \param a, a double
 *  \param A, a SiconosMatrix (IN)
 *  \param B a SiconosMatrix (IN-OUT)
 *  \param init, a bool
 */
void scal(double, const SiconosMatrix&, SiconosMatrix&, bool = true);



#endif
