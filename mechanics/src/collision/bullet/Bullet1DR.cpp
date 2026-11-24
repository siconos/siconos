// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "Bullet1DR.hpp"

#include <BulletCollision/NarrowPhaseCollision/btManifoldPoint.h>

#include "SiconosVector.hpp"

siconos::collision::bullet::Bullet1DR::Bullet1DR(std::shared_ptr<btManifoldPoint> point)
    : siconos::modeling::NewtonEulerNormalVelocityR{}, _contactPoints{point} {}

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
