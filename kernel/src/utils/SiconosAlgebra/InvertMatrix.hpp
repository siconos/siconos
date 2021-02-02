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


/*! \file InvertMatrix.hpp

The following code inverts the matrix input using LU-decomposition with backsubstitution of unit vectors. Reference: Numerical Recipies in C, 2nd ed., by Press, Teukolsky, Vetterling & Flannery.

you can solve Ax=b using three lines of ublas code:

permutation_matrix<> piv;
lu_factorize(A, piv);
lu_substitute(A, piv, x);

*/
#ifndef INVERT_MATRIX_HPP
#define INVERT_MATRIX_HPP

// REMEMBER to update "lu.hpp" header includes from boost-CVS
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

/** Matrix inversion routine.

    Uses lu_factorize and lu_substitute in uBLAS to invert a matrix 
    
    \param input source matrix
    \param output inverted matrix.  

*/
template<class T, class U, class V>
bool InvertMatrix(const boost::numeric::ublas::matrix<T, U, V>& input,
                  boost::numeric::ublas::matrix<T, U, V>& inverse)
{
  typedef boost::numeric::ublas::permutation_matrix<std::size_t> pmatrix;
  // create a working copy of the input
  boost::numeric::ublas::matrix<T, U, V> A(input);
  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());
  
  // perform LU-factorization
  int res = lu_factorize(A,pm);
  if(res != 0) return false;
  
  // create identity matrix of "inverse"
  inverse.assign(boost::numeric::ublas::identity_matrix<T>(A.size1()));
  
// backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);

  return true;
}
#endif //INVERT_MATRIX_HPP


