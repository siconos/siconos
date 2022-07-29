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

/** ! \file SiconosAlgebraTypeDef.hpp
    \brief Header file for Siconos Algebra objects

    This file provides aliases for matrix and vector objects and so on ...
*/

#ifndef SiconosAlgebraTypeDef
#define SiconosAlgebraTypeDef

// Make sure Fortran function have a calling convention compatible with gfortran
#ifndef BIND_FORTRAN_LOWERCASE_UNDERSCORE
#define BIND_FORTRAN_LOWERCASE_UNDERSCORE
#endif

// We do not want to link to any Boost lib, we use header only parts
#ifndef BOOST_ALL_NO_LIB
#define BOOST_ALL_NO_LIB
#endif

//#include <array>
//#include <boost/numeric/ublas/fwd.hpp>
//#include <complex>
// #include <limits>
//#include <vector>

//#include "SiconosConfig.h"

namespace siconos::algebra {

/** siconos::UBLAS_TYPE is an enumerated type of Siconos::DENSE, TRIANGULAR, SYMMETRIC,
    SPARSE, BANDED. It is used to describe the type of matrix or
    vector we want to construct.
*/
enum class UBLAS_TYPE {
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


// /** Complex matrix, column major */
// using complex_matrix =
//     boost::numeric::ublas::matrix<std::complex<double>, boost::numeric::ublas::column_major>;

// Notes:
// Vector definition in boost: vector<T,A> see
// http://www.boost.org/libs/numeric/ublas/doc/vector.htm T: value type A: storage type


// /** Complex vector */
// using complex_vector = boost::numeric::ublas::vector<std::complex<double>>;


/** Some containers for vectors - Used for example in Relation to compute y and r
 * when there are 2 DS*/
//using VectorOfBlockVectors = std::vector<std::shared_ptr<BlockVector>>;

////** Some containers for matrices - Used for example to handle the
///    various jacobian in LagrangianDS. */
/// typedef std::vector<std::shared_ptr<SiconosMatrix>> VectorOfMatrices;

/// /** Some containers for matrices - Used for example to handle the
///    various jacobian in LagrangianDS. */
/// typedef std::vector<std::shared_ptr<SimpleMatrix>> VectorOfSMatrices;

/// /** type of object used to save indices */
/// typedef std::vector<std::size_t> Index;

}  // namespace siconos::algebra

// Set this to use lapack::optimal_workspace where required in lapack routines.
#define USE_OPTIMAL_WORKSPACE

#endif
