/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include "Bullet5DR.hpp"
#include <RigidBodyDS.hpp>
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

#include <boost/math/quaternion.hpp>

static
void copyQuatRot(boost::math::quaternion<double>& from, SiconosVector& to)
{
  to(3) = from.R_component_1();
  to(4) = from.R_component_2();
  to(5) = from.R_component_3();
  to(6) = from.R_component_4();
}

static
void copyQuatRot(const SiconosVector& from, boost::math::quaternion<double>& to)
{
  to = boost::math::quaternion<double>(from(3), from(4), from(5), from(6));
}

static
void copyQuatPos(const boost::math::quaternion<double>& from, SiconosVector& to)
{
  to(0) = from.R_component_2();
  to(1) = from.R_component_3();
  to(2) = from.R_component_4();
}

static
void copyQuatPos(const SiconosVector& from, boost::math::quaternion<double>& to)
{
  to = boost::math::quaternion<double>(0, from(0), from(1), from(2));
}

static
void copyQuatPos(const btVector3& from, boost::math::quaternion<double>& to)
{
  to = boost::math::quaternion<double>(0, from.x(), from.y(), from.z());
}

static void copyBtVector3(const btVector3 &from, SiconosVector& to)
{
  to(0) = from.x();
  to(1) = from.y();
  to(2) = from.z();
}

Bullet5DR::Bullet5DR()
  : Contact5DR()
{
}

void Bullet5DR::updateContactPointsFromManifoldPoint(const btPersistentManifold& manifold,
                                                   const btManifoldPoint& point,
                                                   bool flip, double scaling,
                                                   SP::NewtonEulerDS ds1,
                                                   SP::NewtonEulerDS ds2)
{
  // Get new world positions of contact points and calculate relative
  // to ds1 and ds2

  ::boost::math::quaternion<double> rq1, rq2, posa;
  ::boost::math::quaternion<double> pq1, pq2, posb;
  copyQuatPos(*ds1->q(), pq1);
  copyQuatPos(point.getPositionWorldOnA() / scaling, posa);
  copyQuatRot(*ds1->q(), rq1);
  if (ds2)
  {
    copyQuatPos(*ds2->q(), pq2);
    copyQuatPos(point.getPositionWorldOnB() / scaling, posb);
    copyQuatRot(*ds2->q(), rq2);
  }

  if (flip)
  {
    ::boost::math::quaternion<double> tmp = posa;
    posa = posb;
    posb = tmp;
  }

  SiconosVector va(3), vb(3), vn(3);
  if (flip) {
    copyQuatPos((1.0/rq1) * (posb - pq1) * rq1, va);
    if (ds2)
      copyQuatPos((1.0/rq2) * (posa - pq2) * rq2, vb);
    else {
      // If no body2, position is relative to 0,0,0
      copyBtVector3(point.getPositionWorldOnA() / scaling, vb);
    }
  } else {
    copyQuatPos((1.0/rq1) * (posa - pq1) * rq1, va);
    if (ds2)
      copyQuatPos((1.0/rq2) * (posb - pq2) * rq2, vb);
    else {
      // If no body2, position is relative to 0,0,0
      copyBtVector3(point.getPositionWorldOnB() / scaling, vb);
    }
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

  Contact5DR::updateContactPoints(va, vb, vn*(flip?-1:1));
}
