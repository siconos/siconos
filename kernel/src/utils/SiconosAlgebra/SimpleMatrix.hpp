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

/*! \file SimpleMatrix.hpp
 */

#ifndef __SimpleMatrix__
#define __SimpleMatrix__

// #include <iosfwd>  // for ostream

#include "SiconosMatrix.hpp"         // for SiconosMatrix, MATRIX_UBL...
#include "SiconosVector.hpp"
#include "SiconosSerialization.hpp"  // For ACCEPT_SERIALIZATION

namespace siconos::algebra {

/**
   Matrix (embedded various types of Boost matrices of double)

   SimpleMatrix is used in the platform to store matrices (mathematical object) of double.

   Possible types: Siconos::DENSE (default),
   TRIANGULAR, SYMMETRIC, SPARSE, BANDED, ZERO,
   Siconos::IDENTITY,
   Siconos::SPARSE_COORDINATE.

   \todo: review resize function for Banded, Symetric and Triangular. Error in tests.

*/
using SimpleMatrix = SiconosMatrix;

void normInfByColumn(const SimpleMatrix &m, SiconosVector &v);

void solveInPlace(SimpleMatrix &A, SimpleMatrix &B);

void solveInPlace(SimpleMatrix &A, SimpleMatrix &&B);

}  // namespace siconos::algebra

#endif
