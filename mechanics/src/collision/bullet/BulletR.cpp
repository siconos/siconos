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

#include <BulletCollision/NarrowPhaseCollision/btManifoldPoint.h>
#include <BulletCollision/NarrowPhaseCollision/btPersistentManifold.h>

#include <iostream>

#include "BulletSiconosCommon.hpp"  // for copyQuatPos etc
// #include "siconos_debug.h"

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
}

void siconos::collision::bullet::BulletR::display() const {
  std::cout << "BulletR display()" << std::endl;
  ContactR::display();

  std::cout << "&btObject[0]" << &btObject[0] << std::endl;
  std::cout << "&btObject[1]" << &btObject[1] << std::endl;
  std::cout << "&btShape[0]" << &btShape[0] << std::endl;
  std::cout << "&btShape[1]" << &btShape[1] << std::endl;
}
