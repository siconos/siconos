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

// All the boost bindings required includes ...
#include <boost/numeric/bindings/lapack.hpp>
//#include <boost/numeric/bindings/noop.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

#include "KernelConfig.h"

// Typedef to define interface to boost ublas.
//#include "SiconosAlgebraTypeDef.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

// Some utilities (print ...)
#include "utils.h"

namespace lapack = boost::numeric::bindings::lapack;

namespace Siconos {
  namespace eigenproblems {

    int syev(SiconosVector& eigenval, SiconosMatrix& eigenvec)
    {
      int info = 0;
      // Eigenvec must contains the values of the matrix from which we want
      // to compute eigenvalues and vectors. It must be a symmetric matrix.
      // It will be overwritten with eigenvectors.  
      
      // Adaptor to symmetric_mat. Warning : no copy, eigenvec will be modified
      // by syev.
      symmetric_type s_a(*eigenvec.dense());
#ifdef USE_OPTIMAL_WORKSPACE
      info += lapack::syev('V', s_a, *eigenval.dense(), lapack::optimal_workspace());
#endif
#ifdef USE_MINIMAL_WORKSPACE
      info += lapack::syev('V', s_a, *eigenval.dense(), lapack::minimal_workspace());
#endif
      std::cout << "Compute eigenvalues ..." << std::endl;
      return info;
    }

    int geev(complex_matrix& input_mat, complex_vector& eigenval, complex_matrix& left_eigenvec, complex_matrix& right_eigenvec)
    {
      // input_mat must contains the values of the matrix from which we want
      // to compute eigenvalues and vectors. It must be a complex matrix.
      // It will be overwritten with temp results.  
      int info = 0;
      
#ifdef USE_OPTIMAL_WORKSPACE
      info += lapack::geev('N','V', input_mat, eigenval, left_eigenvec, right_eigenvec, lapack::optimal_workspace());
#endif
#ifdef USE_MINIMAL_WORKSPACE
      info += lapack::geev('N','V', input_mat, eigenval, left_eigenvec, right_eigenvec, lapack::minimal_workspace());
#endif
      return info;
    }

  } // namespace eigenproblems
} // namespace Siconos
