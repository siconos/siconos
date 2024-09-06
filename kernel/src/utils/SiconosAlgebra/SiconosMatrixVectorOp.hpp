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
#include "SiconosMatrix.hpp"
#include <complex>
#include "SiconosVector.hpp"
#include "BlockVector.hpp"

namespace siconos::algebra {


/** compute the product m1 * trans(m2)
 *  \param 2 SiconosVectors
 *  \return a SiconosMatrix
 */
SiconosMatrix outer_prod(const SiconosVector&, const SiconosVector&);

/** prod(A, x, y, init) computes y = A*x or y += A*x if init = false
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y, bool init = true);
void prod(const MapType& A, const SiconosVector& x, SiconosVector& y, bool init = true);

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

/** prod(a, A, x, y, init) computes y = a*A*x or y += a*A*x if init = false
  \param a a double
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param init a bool (default = true)
  */
void prod(double a, const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y,
          bool init = true);

/** subprod(A, x, y) computes sub_y = sub_A*sub_x or sub_y += sub_A*sub_x if init = false
  \param A a SiconosMatrix
  \param x a SiconosVector
  \param[in,out] y a SiconosVector
  \param coord an std::vector<std::size_t> = [r0A r1A c0A c1A r0x r1x r0y r1y];
  subA is the sub-matrix of A, for row numbers between r0A and r1A-1 and columns between c0A
  and c1A-1; The same for x and y with rix and riy. \param init a bool (default = true)
  */
void subprod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y,
             const std::vector<std::size_t>& coord, bool init = true);

void subprod(const SiconosMatrix& A, const BlockVector& x, SiconosVector& y,
             const std::vector<std::size_t>& coord, bool init = true);

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

/** Compute eigenvalues and eigenvectors of a real symmetric matrix A
 *  See examples of use in test/EigenProblemsTest.cpp.
 *  \param[in,out] eigenval : eigenvalues of the matrix
 *  \param[in,out] eigenvec : input matrix A, replace with eigenvectors (columns) in output.
 *  \param[in] withVect : true if eigenvectors are to be computed (default = true).
 *  \return int : return value from lapack routine. 0 if successful.
 */
int syev(SiconosVector& eigenval, SiconosMatrix& eigenvec, bool withVect = true);

/** Compute eigenvalues and eigenvectors of a nonsymmetrix complex matrix
 *  See examples of use in test/EigenProblemsTest.cpp.
 *  \param[in,out] input_mat SiconosMatrix : input matrix.
 *  \param[in,out] eigenval complex_vector : eigenvalues of the matrix
 *  \param[in,out] left_eigenvec complex_matrix : matrix of the left eigenvectors
 *  \param[in, out] right_eigenvec  complex_matrix : matrix of the right eigenvectors
 *  \param[in] withLeft : true if left eigenvectors are to be computed (default = false).
 *  \param[in] withRight : true if right  eigenvectors are to be computed (default = true).
 *  \return int : return value from lapack routine. 0 if succesful.
 */
// int geev(SiconosMatrix& input_mat,
//          boost::numeric::ublas::vector<std::complex<double>>& eigenval,
//          boost::numeric::ublas::matrix<std::complex<double>,
//                                        boost::numeric::ublas::column_major>& left_eigenvec,
//          boost::numeric::ublas::matrix<std::complex<double>,
//                                        boost::numeric::ublas::column_major>& right_eigenvec,
//          bool withLeft = false, bool withRight = true);

}  // namespace siconos::algebra

#endif
