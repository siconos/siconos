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

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES 1
#include "ContactR.hpp"

#include <iostream>

#include "BodyShapeRecord.hpp"
#include "SiconosVector.hpp"
#include "siconos_debug.h"

void siconos::collision::ContactR::computeh(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  DEBUG_BEGIN("siconos::collision::ContactR::computeh(...)\n");

  // Update contact points and distance if necessary
  NewtonEuler3DR::computeh(q1, q2, y);

  y(0) = distance();

  DEBUG_PRINTF("position on A : %g,%g,%g\n", contactPoint1_(0), contactPoint1_(1),
               contactPoint1_(2));
  DEBUG_PRINTF("position on B : %g,%g,%g\n", contactPoint2_(0), contactPoint2_(1),
               contactPoint2_(2));
  DEBUG_PRINTF("normal on B   : %g,%g,%g\n", nc_(0), nc_(1), nc_(2));

  DEBUG_END("siconos::collision::ContactR::computeh(...)\n");
}

void siconos::collision::ContactR::updateContactPoints(
    const siconos::algebra::SiconosVector3& pos1, const siconos::algebra::SiconosVector3& pos2,
    const siconos::algebra::SiconosVector3& normal) {
  // Copy relative positions
  relPc1_ = pos1;
  relPc2_ = pos2;

  // Update normal
  relNc_ = normal;

  assert(!(relNc_(0) == 0 && relNc_(1) == 0 && relNc_(2) == 0) && "nc = 0, problems..\n");
}

void siconos::collision::ContactR::display() const {
  std::cout << "ContactR display()\n";
  if (bodyShapeRecordA) {
    std::cout << "bodyShapeRecordA: \n";
    bodyShapeRecordA->display();
  }
  if (bodyShapeRecordB) {
    std::cout << "bodyShapeRecordB: \n";
    bodyShapeRecordB->display();
  }
}
