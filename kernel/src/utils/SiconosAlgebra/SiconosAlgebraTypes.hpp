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

/** ! \file SiconosAlgebraTypes.hpp
    
    Enum to list matrices and vectors types (dense, sparse ...)
*/

#ifndef SiconosAlgebraTypes
#define SiconosAlgebraTypes

// Make sure Fortran function have a calling convention compatible with gfortran
#ifndef BIND_FORTRAN_LOWERCASE_UNDERSCORE
#define BIND_FORTRAN_LOWERCASE_UNDERSCORE
#endif

// We do not want to link to any Boost lib, we use header only parts
#ifndef BOOST_ALL_NO_LIB
#define BOOST_ALL_NO_LIB
#endif

namespace siconos::algebra {

/** siconos::UblasType is an enumerated type of Siconos::DENSE, TRIANGULAR, SYMMETRIC,
    SPARSE, BANDED. It is used to describe the type of matrix or
    vector we want to construct.
*/
enum class UblasType {
  BLOCK = 0,
  /** id for dense matrix or vector */
  DENSE = 1,
  /** id for triangular matrix */
  TRIANGULAR = 2,
  /** id for symmetric matrix */
  SYMMETRIC = 3,
  /** id for sparse matrix or vector */
  SPARSE = 4,
  /** id for banded matrix */
  BANDED = 5,
  /** id for zero matrix */
  ZERO = 6,
  /** id for identity matrix */
  IDENTITY = 7,
  /** id for sparse matrix or vector */
  SPARSE_COORDINATE = 8
};
}  // namespace siconos::algebra

// Set this to use lapack::optimal_workspace where required in lapack routines.
#define USE_OPTIMAL_WORKSPACE

#endif
