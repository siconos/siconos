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


// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES 1
#include "siconos_debug.h"

#include "BulletR.hpp"
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

#include "BulletSiconosCommon.hpp"

#include "BodyShapeRecord.hpp"

BulletR::BulletR()
  : ContactR()
{
}

void BulletR::updateContactPointsFromManifoldPoint(
  const btPersistentManifold& manifold,
  const btManifoldPoint& point,
  bool flip, double scaling,
  SP::NewtonEulerDS ds1,
  SP::NewtonEulerDS ds2)
{
  if (flip)
  {
    copyBtVector3(-1.0*point.m_normalWorldOnB, *_Nc);
    copyBtVector3(point.getPositionWorldOnA() / scaling, *_Pc2);
    copyBtVector3(point.getPositionWorldOnB() / scaling, *_Pc1);
  }
  else
  {
    copyBtVector3(point.m_normalWorldOnB, *_Nc);
    copyBtVector3(point.getPositionWorldOnA() / scaling, *_Pc1);
    copyBtVector3(point.getPositionWorldOnB() / scaling, *_Pc2);
  }
}

void BulletR::display() const
{
  std::cout << "BulletR display()" << std::endl;
  ContactR::display();


  std::cout << "&btObject[0]" << &btObject[0] << std::endl;
  std::cout << "&btObject[1]" << &btObject[1] << std::endl;
  std::cout << "&btShape[0]" << &btShape[0] << std::endl;
  std::cout << "&btShape[1]" << &btShape[1] << std::endl;
}
