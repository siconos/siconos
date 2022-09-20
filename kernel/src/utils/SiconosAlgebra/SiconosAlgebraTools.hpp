/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

/*! \file SiconosAlgebraTools.hpp
  \brief Standalone functions used with matrices and vectors.
*/
#ifndef SICONOSALGEBRATOOLS_H
#define SICONOSALGEBRATOOLS_H

#include <boost/numeric/ublas/fwd.hpp>
#include <complex>
#include <random>
namespace siconos::algebra {

class BlockVector;
class SiconosVector;
class SiconosMatrix;

/** test if two BlockVectors have the same number of blocks with
    blocks of the same size when at the same position
    \param v1 first vector to compare with
    \param v2 second vecstor to compare with
*/
bool isComparableTo(const BlockVector &v1, const BlockVector &v2);

/** test if two matrices have the same number of blocks with
    blocks of the same dimension when at the same position
    \param v1 first matrix to compare with
    \param v2 second matrix to compare with
*/
bool isComparableTo(const SiconosMatrix &m1, const SiconosMatrix &m2);

/** Compute the matrix exponential Exp = exp(A) for general matrices,
    using scaling and Padé approximation. See expm.hpp.
    \param A : input matrix
    \param Exp : result = exp(A)
    \param computeAndAdd : if true, result = result + exp(A)
**/
void expm(SiconosMatrix &A, SiconosMatrix &Exp, bool computeAndAdd = false);

/** Compute eigenvalues and eigenvectors of a real symmetric matrix A
 *  See examples of use in test/EigenProblemsTest.cpp.
 *  \param[in,out] eigenval : eigenvalues of the matrix
 *  \param[in,out] eigenvec : input matrix A, replace with eigenvectors (columns) in output.
 *  \param[in] withVect : true if eigenvectors are to be computed (default = true).
 *  \return int : return value from lapack routine. 0 if successful.
 */
int syev(SiconosVector &eigenval, SiconosMatrix &eigenvec, bool withVect = true);

// /** Compute eigenvalues and eigenvectors of a nonsymmetrix complex matrix
//  *  See examples of use in test/EigenProblemsTest.cpp.
//  *  \param[in,out] input_mat SiconosMatrix : input matrix.
//  *  \param[in,out] eigenval complex_vector : eigenvalues of the matrix
//  *  \param[in,out] left_eigenvec complex_matrix : matrix of the left eigenvectors
//  *  \param[in, out] right_eigenvec  complex_matrix : matrix of the right eigenvectors
//  *  \param[in] withLeft : true if left eigenvectors are to be computed (default = false).
//  *  \param[in] withRight : true if right  eigenvectors are to be computed (default = true).
//  *  \return int : return value from lapack routine. 0 if succesful.
//  */
int geev(
    SiconosMatrix &input_mat, boost::numeric::ublas::vector<std::complex<double>> &eigenval,
    boost::numeric::ublas::matrix<std::complex<double>, boost::numeric::ublas::column_major>
        &left_eigenvec,
    boost::numeric::ublas::matrix<std::complex<double>, boost::numeric::ublas::column_major>
        &right_eigenvec,
    bool withLeft = false, bool withRight = true);

namespace internal {

template <typename T>
struct RndIntGen {
  RndIntGen(T l, T h) : low(l), high(h) {}

  double operator()() { return dist(gen); }

 private:
  T low{0};
  T high{100};
  std::random_device rd;   // non-deterministic generator
  std::mt19937 gen{rd()};  // to seed mersenne twister.
  std::uniform_real_distribution<T> dist{low, high};
};

/** Random init of a boost ublas matrix
 */
template <typename M, typename T = M::value_type>
void randomize(M &m, T min = 0., T max = 100.)
{
  // using value_type = typename M::value_type;
  for (auto it = m.begin1(); it != m.end1(); ++it)
    generate(it.begin(), it.end(), RndIntGen<T>(min, max));
}

}  // namespace internal
}  // namespace siconos::algebra

#endif
