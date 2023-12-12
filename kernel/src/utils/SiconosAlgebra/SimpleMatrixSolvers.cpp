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

#include "SiconosConfig.h"

//#define BIND_FORTRAN_LOWERCASE_UNDERSCORE

#include "BlockMatrix.hpp"
#include "NumericsToolsNamespace.h"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "cholesky.hpp"
#include <Eigen/Dense>

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

#ifdef DEBUG_MESSAGES
#include <cs.h>
#endif

// namespace lapack = boost::numeric::bindings::lapack;

void siconos::algebra::solveInPlace(SimpleMatrix &A, SiconosVector &B)
{
  if(A.isFactorized == false) {
    A.lu_siconos = new Eigen::FullPivLU<SiconosMatrix>(A);
    A.isFactorized = true;
    B = A.lu_siconos->solve(B); // TODO : avoid temp copy
    // std::cout << "solve results " << std::endl << tmp << std::endl;
  } else {
    B = A.lu_siconos->solve(B);
  }
}

void siconos::algebra::solveInPlace(SimpleMatrix &A, SimpleMatrix &B)
{
  assert(A.rows() == B.rows());
  // std::cout << "A " << std::endl << A << std::endl;

  if(A.isFactorized == true) {
    A.lu_siconos = new Eigen::FullPivLU<SiconosMatrix>(A);
    A.isFactorized = true;
    B = A.lu_siconos->solve(B); // TODO : avoid temp copy
    // std::cout << "solve results " << std::endl << tmp << std::endl;
  } else {
    SimpleMatrix tmp = A.lu().solve(B);
    B = A.lu_siconos->solve(B);
  }
  // std::cout << "A  apres" << std::endl << A << std::endl;
}