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

// TODO: "point" parameter used only by BulletSpaceFilter
BulletR::BulletR(const btManifoldPoint &point,
                 SP::SiconosVector q1, SP::SiconosVector q2,
                 bool flip,
                 double y_correction_A,
                 double y_correction_B,
                 double scaling)
  : ContactR(q1, q2, flip, y_correction_A, y_correction_B, scaling)
{
}

void BulletR::updateContactPoints(const btManifoldPoint& point,
                                  SP::NewtonEulerDS ds1, SP::NewtonEulerDS ds2)
{
  // Get new positions
  btVector3 posa = point.getPositionWorldOnA();
  btVector3 posb = point.getPositionWorldOnB();
  SiconosVector va(3), vb(3), vn(3);
  copyBtVector3(posa, va);
  copyBtVector3(posb, vb);

  // Get new normal
  copyBtVector3(point.m_normalWorldOnB, vn);

  ContactR::updateContactPoints(va, vb, point.getDistance(), vn, ds1, ds2);
}
