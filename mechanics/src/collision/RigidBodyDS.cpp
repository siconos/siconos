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
#include "StorageTools.hpp"

// siconos::collision::RigidBodyDS::RigidBodyDS(
//     Eigen::Ref<siconos::algebra::SiconosVector7> position,
//     Eigen::Ref<siconos::algebra::SiconosVector6> twist, double mass,
//     Eigen::Ref<siconos::algebra::SiconosMatrix33> inertia, siconos::algebra::AliasTag)
//     : siconos::modeling::NewtonEulerDS{position, twist, mass, inertia,
//                                        siconos::algebra::alias_t},
//       _contactors(std::make_shared<siconos::collision::SiconosContactorSet>()) {}

siconos::collision::RigidBodyDS::RigidBodyDS(const siconos::algebra::SiconosVector7& position,
                                             const siconos::algebra::SiconosVector6& twist,
                                             double mass,
                                             const siconos::algebra::SiconosMatrix33& inertia)
    : siconos::modeling::NewtonEulerDS{position, twist, mass, inertia,
                                       siconos::algebra::copy_t},
      _contactors(std::make_shared<siconos::collision::SiconosContactorSet>()) {}

void siconos::collision::RigidBodyDS::compute_extrapolated_position(
    double extrapolationCoefficient) {
  // we compute an extrapolation of the position
  if (!_qExtrapolated)
    _qExtrapolated = std::make_shared<siconos::algebra::SiconosVector>(state_q_->size());

  auto velocityIncrement = std::make_shared<siconos::algebra::SiconosVector>(twist_->size());

  (*_qExtrapolated)(0) = (*velocityIncrement)(0);
  (*_qExtrapolated)(1) = (*velocityIncrement)(1);
  (*_qExtrapolated)(2) = (*velocityIncrement)(2);
  siconos::geometry::quaternionFromTwistVector(*velocityIncrement, *_qExtrapolated);
  const siconos::algebra::SiconosVector7& qold = qMemory().getSiconosVector(0);
  siconos::geometry::compositionLawLieGroup(qold, *_qExtrapolated);
}
