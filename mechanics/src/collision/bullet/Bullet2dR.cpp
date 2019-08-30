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
// #define DEBUG_MESSAGES 
#include <debug.h>

#include "Bullet2dR.hpp"
#include <RigidBody2dDS.hpp>
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

//static
//void copyQuatRot2d(boost::math::quaternion<double>& from, SiconosVector& to)
//{
//  double angle = 2.0 * acos(from.R_component_1());
//  to(2) = angle;
//}

static
void copyQuatRot2d(const SiconosVector& from, boost::math::quaternion<double>& to)
{
  double half_angle = from(2)/2.0;

  to = boost::math::quaternion<double>(cos(half_angle), 0.0, 0.0, sin(half_angle));
}

static
void copyQuatPos2d(const boost::math::quaternion<double>& from, SiconosVector& to)
{
  to(0) = from.R_component_2();
  to(1) = from.R_component_3();
}

static
void copyQuatPos2d(const SiconosVector& from, boost::math::quaternion<double>& to)
{
  to = boost::math::quaternion<double>(0, from(0), from(1), 0.0);
}

static
void copyQuatPos2d(const btVector3& from, boost::math::quaternion<double>& to)
{
  to = boost::math::quaternion<double>(0, from.x(), from.y(), 0.0);
}

static void copyBtVector32d(const btVector3 &from, SiconosVector& to)
{
  to(0) = from.x();
  to(1) = from.y();
}

Bullet2dR::Bullet2dR()
  : Contact2dR()
{
}

void Bullet2dR::updateContactPointsFromManifoldPoint(const btPersistentManifold& manifold,
                                                     const btManifoldPoint& point,
                                                     bool flip, double scaling,
                                                     SP::RigidBody2dDS ds1,
                                                     SP::RigidBody2dDS ds2)
{
  DEBUG_BEGIN("Bullet2dR::updateContactPointsFromManifoldPoint(...)\n");
  // Get new world positions of contact points and calculate relative
  // to ds1 and ds2

  DEBUG_PRINTF("point.getPositionWorldOnA().x() = %8.5e\t", point.getPositionWorldOnA().x());
  DEBUG_PRINTF("point.getPositionWorldOnA().y() = %8.5e\t", point.getPositionWorldOnA().y());
  DEBUG_PRINTF("point.getPositionWorldOnA().z() = %8.5e\n", point.getPositionWorldOnA().z());

  DEBUG_PRINTF("point.getPositionWorldOnB().x() = %8.5e\t", point.getPositionWorldOnB().x());
  DEBUG_PRINTF("point.getPositionWorldOnB().y() = %8.5e\t", point.getPositionWorldOnB().y());
  DEBUG_PRINTF("point.getPositionWorldOnB().z() = %8.5e\n", point.getPositionWorldOnB().z());

  DEBUG_PRINTF("point.m_normalWorldOnB.x() = %8.5e\t", point.m_normalWorldOnB.x());
  DEBUG_PRINTF("point.m_normalWorldOnB.y() = %8.5e\t", point.m_normalWorldOnB.y());
  DEBUG_PRINTF("point.m_normalWorldOnB.z() = %8.5e\n", point.m_normalWorldOnB.z());


  ::boost::math::quaternion<double> rq1, rq2, posa;
  ::boost::math::quaternion<double> pq1, pq2, posb;


  /* we build a unit quaternion from the angle */
  DEBUG_EXPR(ds1->q()->display(););
  copyQuatPos2d(*ds1->q(), pq1);

  copyQuatPos2d(point.getPositionWorldOnA() / scaling, posa);

  copyQuatRot2d(*ds1->q(), rq1);


  if (ds2)
  {
    DEBUG_EXPR(ds2->q()->display(););
    copyQuatPos2d(*ds2->q(), pq2);
    copyQuatPos2d(point.getPositionWorldOnB() / scaling, posb);
    copyQuatRot2d(*ds2->q(), rq2);
  }

  if (flip)
  {
    ::boost::math::quaternion<double> tmp = posa;
    posa = posb;
    posb = tmp;
  }

  SiconosVector va(2), vb(2);
  if (flip)
  {
    copyQuatPos2d((1.0/rq1) * (posb - pq1) * rq1, va);
    if (ds2)
      copyQuatPos2d((1.0/rq2) * (posa - pq2) * rq2, vb);
    else
    {
      // If no body2, position is relative to 0,0,0
      copyBtVector32d(point.getPositionWorldOnA() / scaling, vb);
    }
  }
  else
  {
    copyQuatPos2d((1.0/rq1) * (posa - pq1) * rq1, va);
    if (ds2)
      copyQuatPos2d((1.0/rq2) * (posb - pq2) * rq2, vb);
    else {
      // If no body2, position is relative to 0,0,0
      copyBtVector32d(point.getPositionWorldOnB() / scaling, vb);
    }
  }



  SiconosVector vn(3);
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
    copyBtVector32d(point.m_normalWorldOnB, vn);

  vn.resize(2);
  DEBUG_EXPR(va.display(););
  DEBUG_EXPR(vb.display(););
  DEBUG_EXPR(vn.display(););
  Contact2dR::updateContactPoints(va, vb, vn*(flip?-1:1));
  DEBUG_END("Bullet2dR::updateContactPointsFromManifoldPoint(...)\n");
}
