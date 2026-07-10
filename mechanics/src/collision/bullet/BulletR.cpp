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
#include "BulletR.hpp"

#include <BulletCollision/CollisionDispatch/btCollisionObject.h>
#include <BulletCollision/NarrowPhaseCollision/btManifoldPoint.h>
#include <BulletCollision/NarrowPhaseCollision/btPersistentManifold.h>

#include <boost/math/quaternion.hpp>
#include <iostream>

#include "BulletSiconosCommon.hpp"  // for copyQuatPos etc
#include "NewtonEulerDS.hpp"
#include "RotationQuaternion.hpp"

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_NOCOLOR
#include "siconos_debug.h"

void siconos::collision::bullet::BulletR::updateContactPointsFromManifoldPoint(
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
  DEBUG_PRINTF("position on A : %g,%g,%g\n", contactPoint1_(0), contactPoint1_(1),
               contactPoint1_(2));
  DEBUG_PRINTF("position on B : %g,%g,%g\n", contactPoint2_(0), contactPoint2_(1),
               contactPoint2_(2));
  DEBUG_PRINTF("normal on B   : %g,%g,%g\n", nc_(0), nc_(1), nc_(2));
}

void siconos::collision::bullet::BulletR::display() const {
  std::cout << "BulletR display()" << std::endl;
  ContactR::display();

  std::cout << "&btObject[0]" << &btObject[0] << std::endl;
  std::cout << "&btObject[1]" << &btObject[1] << std::endl;
  std::cout << "&btShape[0]" << &btShape[0] << std::endl;
  std::cout << "&btShape[1]" << &btShape[1] << std::endl;
}

void siconos::collision::bullet::BulletR::updateRelativeContactPointsFromManifoldPoint(
    const btPersistentManifold& manifold, const btManifoldPoint& point, bool flip,
    double scaling, std::shared_ptr<siconos::modeling::NewtonEulerDS> ds1,
    std::shared_ptr<siconos::modeling::NewtonEulerDS> ds2) {
  // Get new world positions of contact points and calculate relative
  // to ds1 and ds2

  ::boost::math::quaternion<double> rq1, rq2, posa;
  ::boost::math::quaternion<double> pq1, pq2, posb;
  siconos::geometry::extractPositionToQuaternion(*ds1->q(), pq1);
  siconos::collision::bullet::copyQuatPos(point.getPositionWorldOnA() / scaling, posa);
  siconos::geometry::extractRotationQuaternion(*ds1->q(), rq1);
  if (ds2) {
    siconos::geometry::extractPositionToQuaternion(*ds2->q(), pq2);
    siconos::collision::bullet::copyQuatPos(point.getPositionWorldOnB() / scaling, posb);
    siconos::geometry::extractRotationQuaternion(*ds2->q(), rq2);
  }

  if (flip) {
    ::boost::math::quaternion<double> tmp = posa;
    posa = posb;
    posb = tmp;
  }

  siconos::algebra::SiconosVector3 va, vb, vn;
  if (flip) {
    siconos::geometry::extractVectorFromQuaternion((1.0 / rq1) * (posb - pq1) * rq1, va);
    if (ds2)
      siconos::geometry::extractVectorFromQuaternion((1.0 / rq2) * (posa - pq2) * rq2, vb);
    else {
      // If no body2, position is relative to 0,0,0
      siconos::collision::bullet::copyBtVector3(point.getPositionWorldOnA() / scaling, vb);
    }
  } else {
    siconos::geometry::extractVectorFromQuaternion((1.0 / rq1) * (posa - pq1) * rq1, va);
    if (ds2)
      siconos::geometry::extractVectorFromQuaternion((1.0 / rq2) * (posb - pq2) * rq2, vb);
    else {
      // If no body2, position is relative to 0,0,0
      siconos::collision::bullet::copyBtVector3(point.getPositionWorldOnB() / scaling, vb);
    }
  }

  // Get new normal
  if (ds2) {
    btQuaternion qn(point.m_normalWorldOnB.x(), point.m_normalWorldOnB.y(),
                    point.m_normalWorldOnB.z(), 0);
    btQuaternion qb1 = manifold.getBody1()->getWorldTransform().getRotation();
    // un-rotate normal into body1 frame
    qn = qb1.inverse() * qn * qb1;
    vn(0) = qn.x();
    vn(1) = qn.y();
    vn(2) = qn.z();
    vn = vn / vn.norm();
  } else
    siconos::collision::bullet::copyBtVector3(point.m_normalWorldOnB, vn);

  ContactR::updateContactPoints(va, vb, vn * (flip ? -1 : 1));
}
