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

#include <Interaction.hpp>

BulletR::BulletR(SP::btManifoldPoint point, bool flip) :
  NewtonEulerFrom3DLocalFrameR(),
  _contactPoints(point),
  _flip(flip)
{
  btVector3 posa = _contactPoints->getPositionWorldOnA();
  btVector3 posb = _contactPoints->getPositionWorldOnB();
  if (flip) {
      posa = _contactPoints->getPositionWorldOnB();
      posb = _contactPoints->getPositionWorldOnA();
  }

  (*pc1())(0) = posa[0];
  (*pc1())(1) = posa[1];
  (*pc1())(2) = posa[2];
  (*pc2())(0) = posb[0];
  (*pc2())(1) = posb[1];
  (*pc2())(2) = posb[2];

  (*nc())(0) = _contactPoints->m_normalWorldOnB[0] * (flip ? -1 : 1);
  (*nc())(1) = _contactPoints->m_normalWorldOnB[1] * (flip ? -1 : 1);
  (*nc())(2) = _contactPoints->m_normalWorldOnB[2] * (flip ? -1 : 1);
}

void BulletR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  DEBUG_BEGIN("BulletR::computeh(...)\n");

  NewtonEulerR::computeh(time, q0, y);

  DEBUG_PRINT("start of computeh\n");

  btVector3 posa = _contactPoints->getPositionWorldOnA();
  btVector3 posb = _contactPoints->getPositionWorldOnB();
  if (_flip) {
      posa = _contactPoints->getPositionWorldOnB();
      posb = _contactPoints->getPositionWorldOnA();
  }

  (*pc1())(0) = posa[0];
  (*pc1())(1) = posa[1];
  (*pc1())(2) = posa[2];
  (*pc2())(0) = posb[0];
  (*pc2())(1) = posb[1];
  (*pc2())(2) = posb[2];

  {
    y.setValue(0, _contactPoints->getDistance());

    (*nc())(0) = _contactPoints->m_normalWorldOnB[0] * (_flip ? -1 : 1);
    (*nc())(1) = _contactPoints->m_normalWorldOnB[1] * (_flip ? -1 : 1);
    (*nc())(2) = _contactPoints->m_normalWorldOnB[2] * (_flip ? -1 : 1);
  }

  DEBUG_PRINTF("distance : %g\n",  y.getValue(0));


  DEBUG_PRINTF("position on A : %g,%g,%g\n", posa[0], posa[1], posa[2]);
  DEBUG_PRINTF("position on B : %g,%g,%g\n", posb[0], posb[1], posb[2]);
  DEBUG_PRINTF("normal on B   : %g,%g,%g\n", (*nc())(0), (*nc())(1), (*nc())(2));

  DEBUG_END("BulletR::computeh(...)\n");


}
