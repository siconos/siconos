/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include <debug.h>

#include "BulletR.hpp"
#include <BodyDS.hpp>
#include <Interaction.hpp>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunreachable-code"
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif

#include <BulletCollision/NarrowPhaseCollision/btManifoldPoint.h>
#include <BulletCollision/CollisionDispatch/btCollisionObject.h>

#include <btBulletCollisionCommon.h>

#if defined(__clang__)
#pragma clang diagnostic pop
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic pop
#endif

static void copyBtVector3(const btVector3 &from, SiconosVector& to)
{
  to(0) = from.x();
  to(1) = from.y();
  to(2) = from.z();
}

BulletR::BulletR()
  : ContactR()
{
}

void BulletR::updateContactPointsFromManifoldPoint(const btPersistentManifold& manifold,
                                                   const btManifoldPoint& point,
                                                   bool flip, double scaling,
                                                   SP::NewtonEulerDS ds1,
                                                   SP::NewtonEulerDS ds2)
{
  // Get new positions
  btVector3 posa = point.m_localPointA * scaling;
  btVector3 posb = point.m_localPointB * scaling;

  SiconosVector va(3), vb(3), vn(3);
  if (flip) {
    if (ds2)
      copyBtVector3(posa, vb);
    else
      // If no body2, position is relative to 0,0,0
      // TODO: scaling?
      copyBtVector3(point.getPositionWorldOnA(), vb);
    copyBtVector3(posb, va);
  } else {
    copyBtVector3(posa, va);
    if (ds2)
      copyBtVector3(posb, vb);
    else
      // If no body2, position is relative to 0,0,0
      // TODO: scaling?
      copyBtVector3(point.getPositionWorldOnB(), vb);
  }

  // Get new normal
  if (ds2)
  {
    btQuaternion qn(point.m_normalWorldOnB.x(),
                    point.m_normalWorldOnB.y(),
                    point.m_normalWorldOnB.z(), 0);
    btQuaternion qb1 = manifold.getBody1()->getWorldTransform().getRotation();
    // un-rotate normal into body1 frame
    qn = qb1.inverse() * qn * qb1;
    vn(0) = qn.x();
    vn(1) = qn.y();
    vn(2) = qn.z();
    vn = vn/vn.norm2();
  }
  else
    copyBtVector3(point.m_normalWorldOnB, vn);

  ContactR::updateContactPoints(va, vb, vn*(flip?-1:1));
}
