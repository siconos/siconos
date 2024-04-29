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
#include "SecondOrderDS.hpp"

#include "BoundaryCondition.hpp"
#include "SiconosVector.hpp"
#include "SiconosMatrix.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <iostream>

#include "siconos_debug.h"

void siconos::modeling::SecondOrderDS::setBoundaryConditions(
    std::shared_ptr<siconos::modeling::BoundaryCondition> newbd) {
  if (_boundaryConditions) {
    std::cout << "Warning : SecondOrderDS::setBoundaryConditions. old boundary conditions "
                 "were pre-existing. \n";
  }
  _boundaryConditions = newbd;
  _reactionToBoundaryConditions =
      std::make_shared<siconos::algebra::SiconosVector>(_boundaryConditions->size());
};

void siconos::modeling::SecondOrderDS::setMassPtr(
    std::shared_ptr<siconos::algebra::SiconosMatrix> newPtr) {
  _mass = newPtr;
  _hasConstantMass = true;
}
