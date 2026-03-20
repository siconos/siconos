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

#include "Bullet1DR.hpp"

#include <BulletCollision/NarrowPhaseCollision/btManifoldPoint.h>

#include "SiconosVector.hpp"

siconos::collision::bullet::Bullet1DR::Bullet1DR(std::shared_ptr<btManifoldPoint> point)
    : siconos::modeling::NewtonEuler1DR{}, _contactPoints{point} {}

void siconos::collision::bullet::Bullet1DR::computeh(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  y(0) = _contactPoints->getDistance();
  auto posa = _contactPoints->getPositionWorldOnA();
  auto posb = _contactPoints->getPositionWorldOnB();
  contactPoint1_ << posa[0], posa[1], posa[2];
  contactPoint1_ << posb[0], posb[1], posb[2];
  nc_ << _contactPoints->m_normalWorldOnB[0], _contactPoints->m_normalWorldOnB[1],
      _contactPoints->m_normalWorldOnB[2];
}
