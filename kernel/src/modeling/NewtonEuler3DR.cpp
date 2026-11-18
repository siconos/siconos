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

#include "NewtonEuler3DR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::modeling::NewtonEuler3DR::initialize(Interaction& inter) {
  const auto qSize = 7 * (inter.getSizeOfDS() / 6);
  H_NE_internal_storage_.resize(3, qSize);
  H_NE_view_ =
      std::make_shared<siconos::algebra::MapType>(H_NE_internal_storage_.data(), 3, qSize);
  H_NE_view_->setZero();
  NewtonEulerR::initialize(inter);
}

void siconos::modeling::NewtonEuler3DR::computeH_NE_prod_T(
    const Interaction& inter, const siconos::algebra::BlockVector& q0) {
  DEBUG_BEGIN("siconos::modeling::NewtonEuler3DR::computeH_NE_prod_T(...)\n");
  computeContactFrame_(true);
  computeContactVelocityJacobian_(q0, 3);
  DEBUG_END("siconos::modeling::NewtonEuler3DR::computeH_NE_prod_T(...)\n");
}
