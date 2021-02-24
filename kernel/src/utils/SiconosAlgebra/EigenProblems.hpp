/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#ifndef EIGENPROBLEMS_HPP
#define EIGENPROBLEMS_HPP
/*! \file EigenProblems.hpp
  Functions to compute eigenvalues/vectors of a matrix.
  Based on boost ublas bindings to lapack.
  Usage : see EigenProblemsTest.cpp.
*/

#include "SiconosMatrix.hpp"

namespace Siconos {
  namespace eigenproblems {

/** Compute eigenvalues and eigenvectors of a real symmetric matrix A
 *   See examples of use in test/EigenProblemsTest.cpp.
 *   \param[in,out] eigenval : eigenvalues of the matrix
 *   \param[in,out] eigenvec : input matrix A, replace with eigenvectors (columns) in output.
 *   \param[in] withVect : true if eigenvectors are to be computed (default = true).
 *   eigenvector. 
 *   \return int : return value from lapack routine. 0 if successful.
 */
    int syev(SiconosVector& eigenval, SiconosMatrix& eigenvec, bool withVect = true);

/** Compute eigenvalues and eigenvectors of a nonsymmetrix complex matrix 
 *   See examples of use in test/EigenProblemsTest.cpp.
 *   \param[in,out] input_mat SiconosMatrix : input matrix.
 *   \param[in,out] eigenval complex_vector : eigenvalues of the matrix
 *   \param[in,out] left_eigenvec complex_matrix : matrix of the left eigenvectors
 *   \param[in, out] right_eigenvec  complex_matrix : matrix of the right eigenvectors
 *   \param[in] withLeft : true if left eigenvectors are to be computed (default = false).
 *   \param[in] withRight : true if right  eigenvectors are to be computed (default = true).
 *    \return int : return value from lapack routine. 0 if succesful.
 */
    int geev(SiconosMatrix& input_mat, complex_vector& eigenval,
             complex_matrix& left_eigenvec, complex_matrix& right_eigenvec,
             bool withLeft = false, bool withRight = true);

  } // namespace eigenproblems
} // namespace Siconos


#endif
