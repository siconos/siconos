/* Siconos-Kernel, Copyright INRIA 2005-2013.
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
