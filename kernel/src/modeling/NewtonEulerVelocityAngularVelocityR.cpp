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

#include "NewtonEulerVelocityAngularVelocityR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "RotationQuaternion.hpp"  // siconos::geometry::computeRotationMatrix
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::modeling::NewtonEulerVelocityAngularVelocityR::initialize(Interaction& inter) {
  const auto qSize = 7 * (inter.getSizeOfDS() / 6);
  constraintJacobian_storage_.resize(5, qSize);
  constraintJacobian_view_ = std::make_shared<siconos::algebra::MapType>(
      constraintJacobian_storage_.data(), 5, qSize);
  constraintJacobian_view_->setZero();
  NewtonEulerR::initialize(inter);
}

void siconos::modeling::NewtonEulerVelocityAngularVelocityR::computeContactMatrix(
    const Interaction& inter, const siconos::algebra::BlockVector& q0) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEulerVelocityAngularVelocityR::computeContactMatrix\n");

  // Rows 0-2: relative velocity at the contact point, in the contact frame
  computeContactFrame_(true);
  computeContactVelocityJacobian_(q0, 3);

  // Rows 3-4: tangential components of the relative angular velocity,
  // [t1 t2]^T (R1 Omega1 - R2 Omega2). The columns of the translational velocities stay zero.
  siconos::algebra::SiconosMatrix33 R;
  for (Eigen::Index body = 0; body < static_cast<Eigen::Index>(q0.numberOfBlocks()); ++body) {
    const double sign = (body == 0) ? 1.0 : -1.0;
    siconos::geometry::computeRotationMatrix(*q0.vector(body), R);
    contactMatrix_->block(3, 6 * body + 3, 2, 3) = sign * contactFrame_.bottomRows(2) * R;
  }

  DEBUG_EXPR(siconos::algebra::print(*contactMatrix_););
  DEBUG_END("siconos::modeling::NewtonEulerVelocityAngularVelocityR::computeContactMatrix\n");
}
