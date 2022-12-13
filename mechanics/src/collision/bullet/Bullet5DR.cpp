/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "Bullet5DR.hpp"

#include <BulletCollision/CollisionDispatch/btCollisionObject.h>
#include <BulletCollision/NarrowPhaseCollision/btManifoldPoint.h>
#include <BulletCollision/NarrowPhaseCollision/btPersistentManifold.h>

#include <boost/math/quaternion.hpp>

#include "BulletSiconosCommon.hpp"  // for copyQuatpos
#include "RigidBodyDS.hpp"
#include "RotationQuaternion.hpp"  // for copyQuatpos
#include "SiconosVector.hpp"

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES 1
#include "siconos_debug.h"

void siconos::collision::bullet::Bullet5DR::updateContactPointsFromManifoldPoint(
    const btPersistentManifold& manifold, const btManifoldPoint& point, bool flip,
    double scaling, std::shared_ptr<siconos::modeling::NewtonEulerDS> ds1,
    std::shared_ptr<siconos::modeling::NewtonEulerDS> ds2)
{
  // Get new world positions of contact points and calculate relative
  // to ds1 and ds2

  boost::math::quaternion<double> rq1, rq2, posa;
  boost::math::quaternion<double> pq1, pq2, posb;
  siconos::geometry::copyQuatPos(*ds1->q(), pq1);
  siconos::geometry::copyQuatRot(*ds1->q(), rq1);
  if (ds2) {
    siconos::geometry::copyQuatPos(*ds2->q(), pq2);
    siconos::geometry::copyQuatRot(*ds2->q(), rq2);
  }

  siconos::collision::bullet::copyQuatPos(point.getPositionWorldOnA() / scaling, posa);
  siconos::collision::bullet::copyQuatPos(point.getPositionWorldOnB() / scaling, posb);

  if (flip) {
    boost::math::quaternion<double> tmp = posa;
    posa = posb;
    posb = tmp;
  }

  siconos::algebra::SiconosVector va{3}, vb{3}, vn{3};
  if (flip) {
    siconos::geometry::copyQuatPos((1.0 / rq1) * (posa - pq1) * rq1, va);
    if (ds2)
      siconos::geometry::copyQuatPos((1.0 / rq2) * (posb - pq2) * rq2, vb);
    else {
      // If no body2, position is relative to 0,0,0
      siconos::collision::bullet::copyBtVector3(point.getPositionWorldOnA() / scaling, vb);
    }
  }
  else {
    siconos::geometry::copyQuatPos((1.0 / rq1) * (posa - pq1) * rq1, va);
    if (ds2)
      siconos::geometry::copyQuatPos((1.0 / rq2) * (posb - pq2) * rq2, vb);
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
    vn = vn / vn.norm2();
  }
  else
    siconos::collision::bullet::copyBtVector3(point.m_normalWorldOnB, vn);

  siconos::collision::Contact5DR::updateContactPoints(va, vb, vn * (flip ? -1 : 1));
}
