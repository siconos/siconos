// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "BulletVelocityAngularVelocityR.hpp"

#include <BulletCollision/NarrowPhaseCollision/btManifoldPoint.h>
#include <BulletCollision/NarrowPhaseCollision/btPersistentManifold.h>

#include "BulletSiconosCommon.hpp"  // for copyQuatpos

void siconos::collision::bullet::BulletVelocityAngularVelocityR::
    updateContactPointsFromManifoldPoint(
        const btPersistentManifold& manifold, const btManifoldPoint& point, bool flip,
        double scaling, std::shared_ptr<siconos::modeling::NewtonEulerDS> ds1,
        std::shared_ptr<siconos::modeling::NewtonEulerDS> ds2) {
  if (flip) {
    siconos::collision::bullet::copyBtVector3(-1.0 * point.m_normalWorldOnB, nc_);
    siconos::collision::bullet::copyBtVector3(point.getPositionWorldOnA() / scaling,
                                              contactPoint2_);
    siconos::collision::bullet::copyBtVector3(point.getPositionWorldOnB() / scaling,
                                              contactPoint1_);
  } else {
    siconos::collision::bullet::copyBtVector3(point.m_normalWorldOnB, nc_);
    siconos::collision::bullet::copyBtVector3(point.getPositionWorldOnA() / scaling,
                                              contactPoint1_);
    siconos::collision::bullet::copyBtVector3(point.getPositionWorldOnB() / scaling,
                                              contactPoint2_);
  }
}
