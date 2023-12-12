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

#include <assert.h>

#include "BlockMatrix.hpp"
#include "SiconosAlgebraTools.hpp"  // for isComparableTo
#include "SiconosAlgebraTypes.hpp"  // UblasType
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"  // For isComparableto
#include "SimpleMatrix.hpp"


void siconos::algebra::scal(double a, const SiconosMatrix &A, SiconosMatrix &B, bool init) {
  // To compute B = a * A (init = true) or B += a*A (init = false).
  if (&A == &B) {
    if (init)
      B *= a;
    else
      B *= (1.0 + a);
  } 
  else {
    if (init)
      B.noalias() = a * A;
    else
      B += a * A;
  }
}
