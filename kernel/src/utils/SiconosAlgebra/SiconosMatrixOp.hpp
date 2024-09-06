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

/*! \file SiconosMatrixOp.hpp
  \brief Toolbox for operators acting on matrices
  For matrix-vector op, see SiconosMatrixVectorOp.hpp
*/
#ifndef SICOMAT_OPH
#define SICOMAT_OPH
#include <memory>
#include "SiconosMatrix.hpp"
#include "BlockMatrix.hpp"
#include "SiconosMatrix.hpp"

namespace siconos::algebra {


// /** multiply the current matrix with a scalar
//  *
//  *  \param m the matrix to operate on
//  *  \param s the scalar
//  *  \return SiconosMatrix&
//  */
// SiconosMatrix& operator*=(SiconosMatrix& m, const double& s);

// SiconosMatrix& operator*=(SiconosMatrix& m, const double& s);

// SiconosMatrix& operator/=(SiconosMatrix& m, const double& s);
// const SiconosMatrix operator/(const SiconosMatrix& A, double a);
// void operator+=(std::shared_ptr<SiconosMatrix> A, std::shared_ptr<SiconosMatrix> B);

// /** multiplication of a matrix by a double
//  *  \param A a SiconosMatrix
//  *  \param a a double
//  *  \return a SiconosMatrix
//  */
// SiconosMatrix operator * (const SiconosMatrix& A, double a);

// /** multiplication of a matrix by a double
//  *  \param a a double
//  *  \param A a SiconosMatrix
//  *  \return a SiconosMatrix
//  */
// SiconosMatrix operator*(double a, const SiconosMatrix& A);

// const SiconosMatrix operator*(const SiconosMatrix&, double);

// // /** multiplication of a SiconosMatrix by a SiconosMatrix
// //  *  \param a const SiconosMatrix&
// //  *  \param a const SiconosMatrix&
// //  *  \return a const SiconosMatrix
// //  */
// // //  const SiconosMatrix operator * (const SiconosMatrix&,const SiconosMatrix&);

// std::shared_ptr<SiconosMatrix> operator*(const std::shared_ptr<SiconosMatrix>,
//                                         const std::shared_ptr<SiconosMatrix>);

// const SiconosMatrix operator/(const SiconosMatrix& A, double a);

// /** multiplication of a std::shared_ptr<SiconosMatrix> by a std::shared_ptr<SiconosMatrix>
//  *  \param A a std::shared_ptr<SiconosMatrix>
//  *  \param B a std::shared_ptr<SiconosMatrix>
//  *  \return a std::shared_ptr<SiconosMatrix>
//  */
// std::shared_ptr<SiconosMatrix> operator * (const std::shared_ptr<SiconosMatrix> A, const
// std::shared_ptr<SiconosMatrix> B);

// /** operator += add B to A
//  *  \param[in,out] A a std::shared_ptr<SiconosMatrix>
//  *  \param B a std::shared_ptr<SiconosMatrix>
//  */
// void operator +=(std::shared_ptr<SiconosMatrix> A, std::shared_ptr<SiconosMatrix> B);

// /** division of the matrix by a double
//  *  \param A a SiconosMatrix
//  *  \param a a double
//  *  \return a SiconosMatrix
//  */
// const SiconosMatrix operator /(const SiconosMatrix& A, double a);

// /** Addition of two matrices, C = A+B
//  * \param A a SiconosMatrix
//  * \param B a SiconosMatrix
//  * \return a SiconosMatrix C
//  */
// const SiconosMatrix operator+(const SiconosMatrix& A, const SiconosMatrix& B);

// /** Addition of two matrices, C = A+B
//  * \param A a std::shared_ptr<SiconosMatrix>
//  * \param B a std::shared_ptr<SiconosMatrix>
//  * \return a std::shared_ptr<SiconosMatrix>
//  */
// std::shared_ptr<SiconosMatrix> operator+(const std::shared_ptr<SiconosMatrix> A,
//                                         const std::shared_ptr<SiconosMatrix> B);

// SiconosMatrix operator+(const SiconosMatrix& A, const SiconosMatrix& B);

/** Addition of two matrices C = A+B
 *  \param A a SiconosMatrix
 *  \param B a SiconosMatrix
 *  \param[in,out] C a SiconosMatrix
 */
void add(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C);

// /** Subtraction of two matrices, C = A-B
//  * \param A a SiconosMatrix
//  * \param B a SiconosMatrix
//  * \return a SiconosMatrix
//  */
// const SiconosMatrix operator-(const SiconosMatrix& A, const SiconosMatrix& B);
// //  const SiconosMatrix operator -(const SiconosMatrix&,const SiconosMatrix&);

/** Subtraction of two matrices C = A-B
 *  \param A a SiconosMatrix
 *  \param B a SiconosMatrix
 *  \param[in,out] C a SiconosMatrix
 */
void sub(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C);

// /**: A==B when (A-B).normInf()<tolerance
//    \param A a SiconosMatrix
//    \param B a SiconosMatrix
//    \return a boolean
//  */
// bool operator==(const SiconosMatrix& A, const SiconosMatrix& B);

// /**: A!=B when (A-B).normInf()>tolerance
//  * \param A a SiconosMatrix
//  * \param B a SiconosMatrix
//  * \return a boolean
//  */
// bool operator!=(const SiconosMatrix& A, const SiconosMatrix& B);

/** multiplication of a matrix by a scalar, B = a*A (init = true) or B += a*A (init = false)
 *  \param a a double
 *  \param A a SiconosMatrix
 *  \param[in,out] B a SiconosMatrix
 *  \param init a bool
 */
void scal(double a, const SiconosMatrix& A, SiconosMatrix& B, bool = true);

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

/** product of two matrices, C = A*B
  \param A a SiconosMatrix
  \param B a SiconosMatrix
  \return C a SiconosMatrix
  */
const SiconosMatrix prod(const SiconosMatrix& A, const SiconosMatrix& B);

/** send data of the matrix to an ostream
 * \param os An output stream
 * \param bm a BlockMatrix
 * \return The same output stream
 */
// std::ostream& operator<<(std::ostream& os, const SiconosMatrix& sm);
// std::ostream& operator<<(std::ostream& os, const SiconosMatrix& sm);
std::ostream& operator<<(std::ostream& os, const BlockMatrix& sm);

// /** Compute the matrix exponential Exp = exp(A) for general matrices,
//   using scaling and Padé approximation. See expm.hpp.
//   \param A : input matrix
//   \param Exp : result = exp(A)
//   \param computeAndAdd : if true, result = result + exp(A)
// **/
// void expm(SiconosMatrix& A, SiconosMatrix& Exp, bool computeAndAdd = false);

/** Copy a subBlock of MIn into a sub-block of MOut - Dim and positions of the sub-block are
 *  given in dim and start.
 *  \param MIn a SPC::SiconosMatrix \param[in,out] MOut a std::shared_ptr<SiconosMatrix>
 *  \param dim an Index, dim[0], dim[1]: number of rows and columns of the sub-block
 *  \param start an Index, start[0], start[1]: position (row, column) of the first
 *  element of the sub-block in MIn start[2], start[3]: position (row, column) of the first
 *  element of the sub-block in MOut.
 */
void setBlock(const SiconosMatrix& MIn, std::shared_ptr<SiconosMatrix> MOut,
              const std::vector<std::size_t>& dim, const std::vector<std::size_t>& start);

/** Copy a subBlock of MIn into a sub-block of MOut - Dim and positions of the sub-block are
 *  given in dim and start.
 *  \param MIn a SPC::SiconosMatrix \param[in,out] MOut a std::shared_ptr<SiconosMatrix>
 *  \param dim an Index, dim[0], dim[1]: number of rows and columns of the sub-block
 *  \param start an Index, start[0], start[1]: position (row, column) of the first
 *  element of the sub-block in MIn start[2], start[3]: position (row, column) of the first
 *  element of the sub-block in MOut.
 */
void setBlock(const SiconosMatrix& MIn, std::shared_ptr<MapType> MOut,
              const std::vector<std::size_t>& dim, const std::vector<std::size_t>& start);


/** test if two matrices have the same number of blocks with
    blocks of the same dimension when at the same position
    \param v1 first matrix to compare with
    \param v2 second matrix to compare with
*/
bool isComparableTo(const SiconosMatrix& m1, const SiconosMatrix& m2);

}  // namespace siconos::algebra

#endif
