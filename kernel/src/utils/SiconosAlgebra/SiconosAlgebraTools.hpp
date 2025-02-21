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

/*! \file SiconosAlgebraTools.hpp
  \brief Standalone functions used with matrices and vectors.
*/
#ifndef SICONOSALGEBRATOOLS_H
#define SICONOSALGEBRATOOLS_H
// <<<<<<< ours
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <random>
namespace siconos::algebra {

/** Matrix inversion routine.
// =======
// #include <algorithm>
// #include <random>
// class SiconosMatrix;
// class BlockVector;
// >>>>>>> theirs

  Uses lu_factorize and lu_substitute in uBLAS to invert a matrix
  Reference: Numerical Recipies in C, 2nd ed., by Press, Teukolsky, Vetterling &
  Flannery.
  \param input source matrix
  \param output inverted matrix.

*/
template <class T, class U, class V>
bool InvertMatrix(const boost::numeric::ublas::matrix<T, U, V> &input,
                  boost::numeric::ublas::matrix<T, U, V> &inverse) {
  // create a working copy of the input
  boost::numeric::ublas::matrix<T, U, V> A(input);
  // create a permutation matrix for the LU-factorization
  boost::numeric::ublas::permutation_matrix<std::size_t> pm(A.size1());

  // perform LU-factorization
  int res = boost::numeric::ublas::lu_factorize(A, pm);
  if (res != 0) return false;

  // create identity matrix of "inverse"
  inverse.assign(boost::numeric::ublas::identity_matrix<T>(A.size1()));

  // backsubstitute to get the inverse
  boost::numeric::ublas::lu_substitute(A, pm, inverse);

  return true;
}
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
template <typename M, typename T = typename M::value_type>
void randomize(M& m, T min = 0., T max = 100.) {
  // using value_type = typename M::value_type;
  for (auto it = m.begin1(); it != m.end1(); ++it)
    std::generate(it.begin(), it.end(), RndIntGen<T>(min, max));
}

}  // namespace internal
}  // namespace siconos::algebra

#endif
