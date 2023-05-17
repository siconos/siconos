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

/*! \file SiconosMatrix.hpp
  Interface for matrices handling.
*/

#ifndef SICOMAT
#define SICOMAT

// #include <boost/numeric/ublas/fwd.hpp>  // boost::numeric fwd
#include "EigenInclude.hpp"
#include <Eigen/Core>
// #include <memory>                       // shared_ptr
// #include <vector>

// #include "CSparseMatrix.h"          // For CSparseMatrix
// #include "SiconosAlgebraTypes.hpp"  // for UblasType
// #include "SiconosException.hpp"
// #include "SiconosSerialization.hpp"  // for ACCEPT_SERIALIZATION

// #include "NumericsFwd.h"  // For NumericsMatrix
// typedef struct NumericsMatrix NumericsMatrix;
struct NumericsMatrix;

namespace siconos::algebra {



/**
   Abstract class to provide interface for matrices handling

   Matrices can be either block or Simple.
   See Derived classes for details.

   In Siconos, a "matrix" can be either a SimpleMatrix or a BlockMatrix, ie a container of
   several pointers to SiconosMatrix

   You can find an overview on how to build and use vectors and matrices in siconos users guide
   .

*/

using SiconosMatrix = Eigen::MatrixXd;

}  // namespace siconos::algebra
#endif
