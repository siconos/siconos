// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "BulletSiconosCommon.hpp"

#include <LinearMath/btVector3.h>

#include <boost/math/quaternion.hpp>

#include "SiconosVector.hpp"

void siconos::collision::bullet::copyQuatPos(const btVector3& from,
                                             boost::math::quaternion<double>& to) {
  to = boost::math::quaternion<double>{0, from.x(), from.y(), from.z()};
}

void siconos::collision::bullet::copyBtVector3(const btVector3& from,
                                               siconos::algebra::SiconosVector3& to) {
  to(0) = from.x();
  to(1) = from.y();
  to(2) = from.z();
}

void siconos::collision::bullet::copyQuatPos2d(const btVector3& from,
                                               boost::math::quaternion<double>& to) {
  to = boost::math::quaternion<double>{0, from.x(), from.y(), 0.0};
}

void siconos::collision::bullet::copyBtVector32d(const btVector3& from,
                                                 siconos::algebra::SiconosVector2& to) {
  to(0) = from.x();
  to(1) = from.y();
}
