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

#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
// All the boost bindings required includes ...
#include <boost/numeric/bindings/lapack.hpp>
// #include <boost/numeric/bindings/noop.hpp>
#include <boost/numeric/bindings/std/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "SiconosConfig.h"
#include "Tools.hpp"

// Typedef to define interface to boost ublas.
// #include "SiconosAlgebraTypeDef.hpp"
#include "EigenProblems.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// Some utilities (print ...)

namespace lapack = boost::numeric::bindings::lapack;

namespace siconos {
namespace eigenproblems {

int syev(SiconosVector& eigenval, SiconosMatrix& eigenvec, bool withVect) {
  int info = 0;
  // Eigenvec must contains the values of the matrix from which we want
  // to compute eigenvalues and vectors. It must be a symmetric matrix.
  // It will be overwritten with eigenvectors.

  // Adaptor to symmetric_mat. Warning : no copy, eigenvec will be modified
  // by syev.

#ifdef USE_OPTIMAL_WORKSPACE
  auto opt = lapack::optimal_workspace();
#endif
#ifdef USE_MINIMAL_WORKSPACE
  auto opt = lapack::minimal_workspace();
#endif

  char jobz;
  if (withVect)
    jobz = 'V';
  else
    jobz = 'N';
  auto num = eigenvec.num();
  if (num == siconos::DENSE) {
    boost::numeric::ublas::symmetric_adaptor<DenseMat, boost::numeric::ublas::lower> s_a(
        *eigenvec.dense());
    info += lapack::syev(jobz, s_a, *eigenval.dense(), opt);
  } else
    THROW_EXCEPTION("Not yet implemented for matrix of type." + enum_to_string(num));

  std::cout << "Compute eigenvalues ..." << std::endl;
  return info;
}

int geev(SiconosMatrix& input_mat, complex_vector& eigenval, complex_matrix& left_eigenvec,
         complex_matrix& right_eigenvec, bool withLeft, bool withRight) {
  int info = 0;
  complex_matrix tmp(*input_mat.dense());
  // tmp must contains the values of the matrix from which we want
  // to compute eigenvalues and vectors. It must be a complex matrix.
  // It will be overwritten with temp results.

  char jobvl, jobvr;
  if (withLeft)
    jobvl = 'V';
  else
    jobvl = 'N';

  if (withRight)
    jobvr = 'V';
  else
    jobvr = 'N';

#ifdef USE_OPTIMAL_WORKSPACE
  info += lapack::geev(jobvl, jobvr, tmp, eigenval, left_eigenvec, right_eigenvec,
                       lapack::optimal_workspace());
#endif
#ifdef USE_MINIMAL_WORKSPACE
  info += lapack::geev(jobvl, jobvr, tmp, eigenval, left_eigenvec, right_eigenvec,
                       lapack::minimal_workspace());
#endif
  return info;
}

}  // namespace eigenproblems
}  // namespace siconos
