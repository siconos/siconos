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
#include "RigidBodyDS.hpp"

#include "RotationQuaternion.hpp"  // for quaternionFromTwistVector and compositionLawLieGroup
#include "SiconosContactor.hpp"
#include "SiconosVector.hpp"
#include "SiconosVisitor.hpp"

siconos::collision::RigidBodyDS::RigidBodyDS(
    Eigen::Ref<siconos::algebra::SiconosVector> position,
    Eigen::Ref<siconos::algebra::SiconosVector> velocity, double mass,
    Eigen::Ref<siconos::algebra::SiconosMatrix> inertia)
    : siconos::modeling::NewtonEulerDS{position, velocity, mass, inertia},
      _contactors(std::make_shared<siconos::collision::SiconosContactorSet>()) {}

void siconos::collision::RigidBodyDS::compute_extrapolated_position(
    double extrapolationCoefficient) {
  // we compute an extrapolation of the position
  if (!_qExtrapolated)
    _qExtrapolated = std::make_shared<siconos::algebra::SiconosVector>(state_q_->size());

  auto velocityIncrement = std::make_shared<siconos::algebra::SiconosVector>(twist_->size());

  _qExtrapolated->setValue(0, velocityIncrement->getValue(0));
  _qExtrapolated->setValue(1, velocityIncrement->getValue(1));
  _qExtrapolated->setValue(2, velocityIncrement->getValue(2));
  siconos::geometry::quaternionFromTwistVector(*velocityIncrement, *_qExtrapolated);
  const auto &qold = qMemory().getSiconosVector(0);
  siconos::geometry::compositionLawLieGroup(qold, *_qExtrapolated);
}
